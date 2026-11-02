"""Second water model (run D) and thermal high-pressure phases (run E).

Run D, TIP4P/2005 versus TIP4P/Ice on the same starting structures:
  D_tip4p2005_<phase>_p<k>_220   Ih, Ic, Isd, liquid at 220 K (same temperature)
  D_tip4p2005_<phase>_p<k>_200   Ih, Ic, liquid at 200 K (same undercooling:
                                 melting points about 252 K and 272 K)
  D_tip4pice_<phase>_p<k>_220    TIP4P/Ice control with the run D protocol,
                                 so that model and protocol effects separate
  Structures are the run C preparations p1-p3 (md_replicas.py).

Run E, TIP4P/Ice high-pressure phases at 220 K, two GenIce preparations each
(different GenIce seed, so different proton configurations for disordered phases):
  E_<phase>_p<k>_220_<P>GPa      II and III at 0.3, V at 0.5, VI at 1.0, VII and
                                 VIII at 3.0 GPa; liquid at 0.3, 0.5, 1.0 and 3.0 GPa
  Ice II and VIII are proton-ordered; the rest proton-disordered. Bulk crystals
  without interfaces stay metastable on these time scales even where another
  phase is stable (e.g. VII rather than VIII at 220 K); retention is checked.

Protocol (differs from runs A-C, which are NVT in the constructed box): minimize;
NPT at the target pressure (crystals 200 ps, Monte Carlo barostat, anisotropic
for orthorhombic and flexible for triclinic cells; liquids are first melted at
300 K for 100 ps, then 1 ns at the target temperature, 2 ns at 200 K); then NVT
in the final NPT box, 50 ps, and 60 frames 5 ps apart, each also FIRE-quenched
for 300 steps. Langevin 1/ps, 2 fs. Per frame we store the potential energy, the
mean oxygen displacement from the first production frame and the mean cosine
between each molecule's dipole direction now and in the first production frame
(phase and proton-order retention checks that use no benchmarked parameter).

Needs OpenMM with CUDA (e.g. the gencompas conda environment) and GenIce2
(.genice-env) for run E:
  python md_conditions.py --list
  python md_conditions.py all | D | E | <runs>
Outputs go to cache/conditions/ (not committed): <run>.npz, <run>.json and the
starting structures.
"""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

import md_replicas

ROOT = Path(__file__).resolve().parents[1]
SEEDTEST = md_replicas.SEEDTEST
OUT = Path(__file__).resolve().parent/"cache/conditions"
GENICE = ROOT/".genice-env/bin/genice2"
sys.path.insert(0, str(SEEDTEST))

FRAMES, STRIDE_PS, NVT_PS = 60, 5, 50
NPT_PS = {"crystal": 200, "liquid": 1000, "liquid_200": 2000}
MELT_PS, MELT_T = 100, 300.
MIN_HEIGHT = 24.  # A; twice the 9 A cutoff plus room for compression
HP_PHASES = {"II": ("2", 0.3), "III": ("3", 0.3), "V": ("5", 0.5), "VI": ("6", 1.0),
             "VII": ("7", 3.0), "VIII": ("8", 3.0)}
HP_LIQUIDS = (0.3, 0.5, 1.0, 3.0)
D_PREPS, E_PREPS = (1, 2, 3), (1, 2)

# Charges (e), LJ sigma (nm) and epsilon (kJ/mol) on O, O-M distance (nm).
MODELS = {
    "tip4pice": dict(q_h=.5897, sigma=.31668, epsilon=.88211, d_om=.01577),
    "tip4p2005": dict(q_h=.5564, sigma=.31589, epsilon=.774907, d_om=.01546),
}
D_OH, HOH = .09572, 104.52


def runs():
    out = [f"D_tip4p2005_{p}_p{k}_220" for p in ("Ih", "Ic", "Isd", "liquid") for k in D_PREPS]
    out += [f"D_tip4p2005_{p}_p{k}_200" for p in ("Ih", "Ic", "liquid") for k in D_PREPS]
    out += [f"D_tip4pice_{p}_p{k}_220" for p in ("Ih", "Ic", "Isd", "liquid") for k in D_PREPS]
    out += [f"E_{p}_p{k}_220_{P}GPa" for p, (_, P) in HP_PHASES.items() for k in E_PREPS]
    out += [f"E_liquid_p{k}_220_{P}GPa" for P in HP_LIQUIDS for k in E_PREPS]
    return out


def parse(run):
    parts = run.split("_")
    if parts[0] == "D":
        _, model, phase, prep, temp = parts
        pressure = 1e-4
    else:
        _, phase, prep, temp, pressure = parts
        model, pressure = "tip4pice", float(pressure.removesuffix("GPa"))
    return dict(dataset=parts[0], model=model, phase=phase, prep=int(prep[1:]), temperature=float(temp),
                pressure_gpa=pressure)


def cell_from_gro(line):
    v = [float(x)*10 for x in line.split()]
    if len(v) == 3:
        return np.diag(v)
    return np.array([[v[0], v[3], v[4]], [v[5], v[1], v[6]], [v[7], v[8], v[2]]])


def reduce_cell(cell):
    """OpenMM/GROMACS reduced form: a along x, b in xy, |b_x| <= a_x/2, |c_x| <= a_x/2, |c_y| <= b_y/2."""
    a, b, c = cell.copy()
    c -= b*round(c[1]/b[1])
    c -= a*round(c[0]/a[0])
    b -= a*round(b[0]/a[0])
    return np.array([a, b, c])


def minimum_image_cell(d, cell):
    frac = d @ np.linalg.inv(cell)
    return (frac - np.round(frac)) @ cell


def heights(cell):
    volume = abs(np.linalg.det(cell))
    return volume/np.linalg.norm(np.cross(cell[[1, 2, 0]], cell[[2, 0, 1]]), axis=1)


def write_gro(path, xyz, cell, title):
    lines = [title, str(4*len(xyz))]
    k = 0
    for m, w in enumerate(xyz):
        for name, p in zip(("OW", "HW1", "HW2", "MW"), w/10):
            k += 1
            lines.append(f"{(m+1) % 100000:5d}{'SOL':<5s}{name:>5s}{k % 100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    c = cell/10
    lines.append(" ".join(f"{x:.5f}" for x in (c[0, 0], c[1, 1], c[2, 2], c[0, 1], c[0, 2], c[1, 0],
                                                  c[1, 2], c[2, 0], c[2, 1])))
    path.write_text("\n".join(lines) + "\n")
    path.with_suffix(".top").write_text(
        '#include "./charmm36-mar2019.ff/forcefield.itp"\n#include "tip4p-sol.itp"\n\n'
        f"[ system ]\n{title}\n\n[ molecules ]\nSOL {len(xyz)}\n")


def genice(lattice, seed, rep):
    command = [str(GENICE), lattice, "--rep", *map(str, rep), "--seed", str(seed),
               "--water", "tip3p", "--format", "gromacs", "--quiet"]
    text = subprocess.run(command, text=True, capture_output=True, check=True,
                          env=dict(os.environ, PYTHONHASHSEED="0")).stdout
    lines = text.splitlines()
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in lines[2:-1]])*10
    return xyz.reshape(-1, 3, 3), cell_from_gro(lines[-1]), command


def prepare_hp(phase, prep):
    """GenIce crystal with every periodic height >= MIN_HEIGHT, as TIP4P (M added)."""
    name = f"E_{phase}_p{prep}"
    gro, meta_path = OUT/f"{name}.gro", OUT/f"{name}.prep.json"
    if gro.exists() and meta_path.exists():
        return gro, json.loads(meta_path.read_text())
    lattice = HP_PHASES[phase][0]
    _, unit_cell, _ = genice(lattice, prep, (1, 1, 1))
    rep = [math.ceil(MIN_HEIGHT/h) for h in heights(unit_cell)]
    waters, cell, command = genice(lattice, prep, rep)
    cell = reduce_cell(cell)
    o, h1, h2 = waters[:, 0], waters[:, 1], waters[:, 2]
    a = MODELS["tip4pice"]["d_om"]/(2*D_OH*np.cos(np.deg2rad(HOH/2)))
    m = o + a*(h1-o) + a*(h2-o)
    write_gro(gro, np.stack([o, h1, h2, m], axis=1), cell, f"ice {phase} GenIce seed {prep}")
    meta = dict(source="genice2", lattice=lattice, genice_seed=prep, rep=rep, command=["genice2", *command[1:]],
                molecules=len(o), min_height=float(heights(cell).min()))
    meta_path.write_text(json.dumps(meta, indent=1))
    return gro, meta


def prepare(info):
    if info["dataset"] == "E" and info["phase"] != "liquid":
        return prepare_hp(info["phase"], info["prep"])
    return md_replicas.prepare(info["phase"], info["prep"])  # run C structures (liquid: seed-test liquid)


def set_model(system, model):
    """Turn the TIP4P/Ice system of system_ice into the requested TIP4P model (O, H, H, M order)."""
    import openmm as mm
    p = MODELS[model]
    a = p["d_om"]/(2*D_OH*np.cos(np.deg2rad(HOH/2)))
    nonbonded = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    for i in range(system.getNumParticles()):
        site = i % 4
        if site == 0:
            nonbonded.setParticleParameters(i, 0., p["sigma"], p["epsilon"])
        elif site in (1, 2):
            nonbonded.setParticleParameters(i, p["q_h"], 1., 0.)
        else:
            nonbonded.setParticleParameters(i, -2*p["q_h"], 1., 0.)
            system.setVirtualSite(i, mm.ThreeParticleAverageSite(i-3, i-2, i-1, 1-2*a, a, a))
    for k in range(nonbonded.getNumExceptions()):  # intramolecular exclusions stay zero
        i, j, *_ = nonbonded.getExceptionParameters(k)
        nonbonded.setExceptionParameters(k, i, j, 0., 1., 0.)
    return system


def dipoles(x):
    o = x[::4]
    d = (x[1::4]-o) + (x[2::4]-o)
    return d/np.linalg.norm(d, axis=1, keepdims=True)


def simulate(run):
    import openmm as mm
    from openmm import app, unit
    import fire
    import system_ice

    info = parse(run)
    phase, temp, model = info["phase"], info["temperature"], info["model"]
    gro_path, meta = prepare(info)
    seed = 1000*info["prep"] + int(temp) + sum(map(ord, run))
    gro = app.GromacsGroFile(str(gro_path))
    top = app.GromacsTopFile(str(gro_path.with_suffix(".top")), periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir=str(SEEDTEST))
    cell0 = np.array(gro.getPeriodicBoxVectors().value_in_unit(unit.angstrom))
    orthorhombic = np.allclose(cell0, np.diag(np.diag(cell0)))
    liquid = phase == "liquid"
    bar = info["pressure_gpa"]*1e4
    platform = mm.Platform.getPlatformByName("CUDA")
    start = time.time()

    def context(t, barostat):
        system, integrator = system_ice.make_system(top, temp=t)
        set_model(system, model)
        if barostat:
            if liquid:
                system.addForce(mm.MonteCarloBarostat(bar*unit.bar, t*unit.kelvin, 25))
            elif orthorhombic:
                system.addForce(mm.MonteCarloAnisotropicBarostat(mm.Vec3(bar, bar, bar)*unit.bar, t*unit.kelvin,
                                                                 True, True, True, 25))
            else:
                system.addForce(mm.MonteCarloFlexibleBarostat(bar*unit.bar, t*unit.kelvin, 25))
        integrator.setRandomNumberSeed(seed)
        ctx = mm.Context(system, integrator, platform, {"Precision": "mixed"})
        return ctx, integrator, system

    def carry(src, dst):
        s = src.getState(getPositions=True)
        dst.setPeriodicBoxVectors(*s.getPeriodicBoxVectors())
        dst.setPositions(s.getPositions())

    def npt(ctx, integrator, ps):
        volumes = []
        for _ in range(ps):
            integrator.step(500)
            volumes.append(ctx.getState().getPeriodicBoxVolume().value_in_unit(unit.angstrom**3))
        return volumes

    ctx, integrator, _ = context(temp, True)
    ctx.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    ctx.setPositions(gro.getPositions())
    ctx.computeVirtualSites()
    mm.LocalEnergyMinimizer.minimize(ctx, 10., 2000)
    if liquid:  # melt with new velocities at the target pressure, then return to temp
        hot, hot_integrator, _ = context(MELT_T, True)
        carry(ctx, hot)
        hot.setVelocitiesToTemperature(MELT_T*unit.kelvin, seed)
        npt(hot, hot_integrator, MELT_PS)
        carry(hot, ctx)
        del hot
    ctx.setVelocitiesToTemperature(temp*unit.kelvin, seed)
    npt_ps = NPT_PS["liquid_200" if liquid and temp < 210 else "liquid" if liquid else "crystal"]
    volumes = npt(ctx, integrator, npt_ps)
    state = ctx.getState(getPositions=True, getVelocities=True)
    del ctx

    nvt, integrator, _ = context(temp, False)
    nvt.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    nvt.setPositions(state.getPositions())
    nvt.setVelocities(state.getVelocities())
    cell = np.array(state.getPeriodicBoxVectors().value_in_unit(unit.angstrom))
    integrator.step(NVT_PS*500)
    frames, energies = [], []
    for _ in range(FRAMES):
        integrator.step(STRIDE_PS*500)
        s = nvt.getState(getPositions=True, getEnergy=True)
        frames.append(s.getPositions(asNumpy=True).value_in_unit(unit.angstrom).astype(np.float32))
        energies.append(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    md_seconds = time.time() - start
    del nvt
    o0, u0 = frames[0][::4].astype(float), dipoles(frames[0].astype(float))
    displacement = [float(np.linalg.norm(minimum_image_cell(x[::4]-o0, cell), axis=1).mean()) for x in frames]
    orientation = [float(np.mean(np.sum(dipoles(x.astype(float))*u0, axis=1))) for x in frames]

    systemq, _ = system_ice.make_system(top)
    set_model(systemq, model)
    fq = fire.fire_integrator()
    cq = mm.Context(systemq, fq, platform, {"Precision": "mixed"})
    cq.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    fire.set_inverse_masses(fq, systemq)
    quenched = []
    for x in frames:
        cq.setPositions(x.astype(float)/10.)
        cq.setVelocitiesToTemperature(0.)
        fq.setGlobalVariableByName("dtc", .001)
        fq.setGlobalVariableByName("alpha", .1)
        fq.setGlobalVariableByName("npos", 0)
        fq.step(300)
        quenched.append(cq.getState(getPositions=True).getPositions(asNumpy=True)
                        .value_in_unit(unit.angstrom).astype(np.float32))
    np.savez_compressed(OUT/f"{run}.npz", x=np.array(frames), xq=np.array(quenched), cell=cell,
                        energy=np.array(energies), displacement=np.array(displacement),
                        orientation=np.array(orientation), npt_volume=np.array(volumes))
    n = len(frames[0])//4
    simulated_ps = npt_ps + NVT_PS + FRAMES*STRIDE_PS + (MELT_PS if liquid else 0)
    out = dict(run=run, **info, structure=gro_path.name, velocity_seed=seed, molecules=n,
               npt_ps=npt_ps, nvt_ps=NVT_PS, frames=FRAMES, stride_ps=STRIDE_PS, simulated_ps=simulated_ps,
               volume_start=volumes[0], volume_end=float(np.mean(volumes[-npt_ps//4:])),
               volume_drift_last_half=float(np.mean(volumes[-npt_ps//4:])/np.mean(volumes[npt_ps//2:3*npt_ps//4])-1),
               density_g_cm3=float(n*18.01528/0.6022140857/abs(np.linalg.det(cell))),
               md_seconds=md_seconds, ns_per_day=simulated_ps/1000/(md_seconds/86400),
               openmm=mm.__version__, **{f"structure_{k}": v for k, v in meta.items()})
    (OUT/f"{run}.json").write_text(json.dumps(out, indent=1))
    print(f"{run}: {md_seconds:.0f} s, density {out['density_g_cm3']:.3f}, "
          f"energy {np.mean(energies[:10])/n:.2f} -> {np.mean(energies[-10:])/n:.2f} kJ/mol/molecule, "
          f"O displacement {displacement[-1]:.2f} A, dipole cosine {orientation[-1]:.2f}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="*")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--force", action="store_true", help="rerun runs that already have output")
    args = parser.parse_args()
    if args.list:
        print("\n".join(runs()))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    md_replicas.OUT.mkdir(parents=True, exist_ok=True)
    selected = runs() if args.runs == ["all"] else [r for r in runs() if r.startswith(args.runs[0] + "_")] \
        if args.runs in (["D"], ["E"]) else args.runs
    for run in selected:
        if (OUT/f"{run}.npz").exists() and not args.force:
            print(f"{run}: exists, skipped")
            continue
        simulate(run)


if __name__ == "__main__":
    main()
