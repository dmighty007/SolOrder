"""Independent TIP4P/Ice preparations for the held-out benchmark (run C).

Each preparation is a new starting structure, not only new velocities:
  Ih, Ic   the same oxygen lattice as the seed test, with a new proton
           configuration (ice-rule loop moves from an Euler-circuit start,
           keeping the net dipole small), run at 220 K and 250 K
  Isd      a new random stacking sequence (5-7 of 12 layers hexagonal) and a
           new proton configuration, 220 K
  liquid   the 220 K liquid of the seed test, melted at 300 K for 50 ps with a
           new velocity seed, then re-equilibrated at 220 K

Protocol (as md_train.py in ice_nucleation_pg, so runs A, B and C compare):
NVT, Langevin 1/ps, 2 fs, constant box; minimize, equilibrate, then 60 frames
5 ps apart, each also FIRE-quenched for 300 steps. Per frame we also store the
potential energy and, for crystals, the mean oxygen displacement from the
minimized lattice (minimum image), as phase-retention checks that use none of
the benchmarked order parameters.

Needs OpenMM with CUDA (e.g. the gencompas conda environment):
  python md_replicas.py --list
  python md_replicas.py Ih_p1_220 Ic_p1_250 ...    # or "all"
Outputs go to cache/replicas/ (not committed): <run>.npz, <run>.gro/.top,
<run>.json (seeds, stacking, dipole, timings).
"""

import argparse
import json
from pathlib import Path
import sys
import time

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

SEEDTEST = Path("/home/dm/Dibyendu/Projects/ice_nucleation_pg/data/seedtest")
OUT = Path(__file__).resolve().parent/"cache/replicas"
sys.path.insert(0, str(SEEDTEST))

# Lattice of ice_nucleation_pg/data/seedtest/build.py (same box and spacing as runs A and B).
B = 2.76
A = B*np.sqrt(8/3)
NX, NY, NL = 8, 5, 12
SHIFT = {0: 0., 1: A/np.sqrt(3), 2: 2*A/np.sqrt(3)}
PREPS = (1, 2, 3, 4, 5)
FRAMES, STRIDE_PS = 60, 5
EQUIL_PS = {"crystal": 50, "liquid": 1000}
MELT_PS, MELT_T = 50, 300.


def runs():
    out = [f"{p}_p{k}_{t}" for p in ("Ih", "Ic") for k in PREPS for t in (220, 250)]
    out += [f"Isd_p{k}_220" for k in PREPS] + [f"liquid_p{k}_220" for k in PREPS]
    return out


def stacking(kind, rng):
    if kind == "Ih":
        return [0, 1]*(NL//2)
    if kind == "Ic":
        return [0, 1, 2]*(NL//3)
    while True:  # random, consecutive layers differ (including the periodic wrap), 5-7 hexagonal
        s = [0]
        for _ in range(NL-1):
            s.append(rng.choice([x for x in range(3) if x != s[-1]]))
        hexagonal = [s[k-1] == s[(k+1) % NL] for k in range(NL)]
        if s[-1] != s[0] and 5 <= sum(hexagonal) <= 7:
            return s


def crystal(sequence):
    box = np.array([NX*A, NY*A*np.sqrt(3), NL*4*B/3])
    pos, layer = [], []
    for k in range(NL):
        for ix in range(NX):
            for iy in range(NY):
                for bx, by in ((0., 0.), (A/2, A*np.sqrt(3)/2)):
                    for dz in (0., B):
                        pos.append((ix*A + bx, iy*A*np.sqrt(3) + by + SHIFT[sequence[k]], k*4*B/3 + dz))
                        layer.append(k)
    hexagonal = np.array([sequence[k-1] == sequence[(k+1) % NL] for k in range(NL)])
    return np.mod(np.array(pos), box), box, np.where(hexagonal[np.array(layer)], 1, 2).astype(np.int8)


def minimum_image(d, box):
    return d - box*np.round(d/box)


def proton_configuration(o, box, rng, sweeps=20, dipole_tol=.02):
    """Donor lists obeying the ice rules, randomized by loop reversals.

    Start: orient every O-O edge along an Euler circuit (two donations per O).
    Move: walk along donated bonds from a random O until the walk revisits an
    oxygen and reverse that loop (keeps 2-in/2-out). A move is accepted unless
    it raises the net dipole above dipole_tol (per molecule, in units of the
    O-O unit vector sum) and above its current value.
    """
    pairs = cKDTree(o, boxsize=box).query_pairs(B*1.1, output_type="ndarray")
    assert np.all(np.bincount(pairs.ravel(), minlength=len(o)) == 4)
    g = nx.MultiGraph()
    g.add_edges_from(map(tuple, pairs))
    out = [[] for _ in o]
    for component in nx.connected_components(g):
        for a, b in nx.eulerian_circuit(g.subgraph(component), source=int(rng.choice(list(component)))):
            out[a].append(b)
    unit = {}
    for a, b in pairs:
        d = minimum_image(o[b]-o[a], box)
        unit[a, b], unit[b, a] = d/np.linalg.norm(d), -d/np.linalg.norm(d)

    def dipole():
        return sum(unit[i, j] for i in range(len(o)) for j in out[i])/len(o)

    m = dipole()
    accepted = 0
    for _ in range(sweeps*len(o)):
        path, seen = [int(rng.integers(len(o)))], {}
        while path[-1] not in seen:
            seen[path[-1]] = len(path)-1
            path.append(out[path[-1]][rng.integers(2)])
        loop = path[seen[path[-1]]:]
        change = sum(unit[b, a]-unit[a, b] for a, b in zip(loop[:-1], loop[1:]))
        new = m + change/len(o)
        if np.linalg.norm(new) > dipole_tol and np.linalg.norm(new) > np.linalg.norm(m):
            continue
        for a, b in zip(loop[:-1], loop[1:]):
            out[a].remove(b)
            out[b].append(a)
        m, accepted = new, accepted+1
    assert all(len(d) == 2 for d in out)
    return out, float(np.linalg.norm(dipole())), accepted


def water_sites(o, box, donors):
    """O, H, H, M positions (A): H at 0.9572 A, HOH 104.52 deg in the donor-bond plane."""
    half = np.deg2rad(104.52)/2
    xyz = np.zeros((len(o), 4, 3))
    for i, (j, k) in enumerate(donors):
        u = [minimum_image(o[x]-o[i], box) for x in (j, k)]
        u = [v/np.linalg.norm(v) for v in u]
        bisector = (u[0]+u[1])/np.linalg.norm(u[0]+u[1])
        normal = (u[0]-u[1])/np.linalg.norm(u[0]-u[1])
        h1 = o[i] + .9572*(np.cos(half)*bisector + np.sin(half)*normal)
        h2 = o[i] + .9572*(np.cos(half)*bisector - np.sin(half)*normal)
        xyz[i] = o[i], h1, h2, o[i] + .13458*(h1-o[i]) + .13458*(h2-o[i])
    return xyz


def write_gro(path, xyz, box, title):
    lines = [title, str(4*len(xyz))]
    k = 0
    for m, w in enumerate(xyz):
        for name, p in zip(("OW", "HW1", "HW2", "MW"), w/10):
            k += 1
            lines.append(f"{(m+1) % 100000:5d}{'SOL':<5s}{name:>5s}{k % 100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    lines.append(f"{box[0]/10:10.5f}{box[1]/10:10.5f}{box[2]/10:10.5f}")
    path.write_text("\n".join(lines) + "\n")
    path.with_suffix(".top").write_text(
        '#include "./charmm36-mar2019.ff/forcefield.itp"\n#include "tip4p-sol.itp"\n\n'
        f"[ system ]\n{title}\n\n[ molecules ]\nSOL {len(xyz)}\n")


def prepare(phase, prep):
    """Starting structure for one preparation; identical for both temperatures."""
    name = f"{phase}_p{prep}"
    gro = OUT/f"{name}.gro"
    meta_path = OUT/f"{name}.prep.json"
    if gro.exists() and meta_path.exists():
        return gro, json.loads(meta_path.read_text())
    if phase == "liquid":
        gro.write_text((SEEDTEST/"liquid.gro").read_text())
        gro.with_suffix(".top").write_text((SEEDTEST/"liquid.top").read_text())
        meta = dict(source="seedtest liquid.gro")
    else:
        rng = np.random.default_rng([2026, sum(map(ord, phase)), prep])
        sequence = stacking(phase, rng)
        o, box, stack = crystal(sequence)
        donors, dipole, accepted = proton_configuration(o, box, rng)
        seq = "".join("ABC"[x] for x in sequence)
        write_gro(gro, water_sites(o, box, donors), box, f"bulk {phase} {seq} prep {prep}")
        np.save(OUT/f"{name}.stack.npy", stack)
        meta = dict(stacking=seq, hexagonal_fraction=float(np.mean(stack == 1)), dipole=dipole,
                    loop_moves=accepted)
    meta_path.write_text(json.dumps(meta, indent=1))
    return gro, meta


def simulate(run):
    import openmm as mm
    from openmm import app, unit
    import fire
    import system_ice

    phase, prep, temp = run.split("_")
    prep, temp = int(prep[1:]), float(temp)
    gro_path, meta = prepare(phase, prep)
    seed = 1000*prep + int(temp) + (7 if phase == "liquid" else sum(map(ord, phase)))
    gro = app.GromacsGroFile(str(gro_path))
    top = app.GromacsTopFile(str(gro_path.with_suffix(".top")), periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir=str(SEEDTEST))
    platform = mm.Platform.getPlatformByName("CUDA")
    start = time.time()

    def context(t):
        system, integrator = system_ice.make_system(top, temp=t)
        integrator.setRandomNumberSeed(seed)
        ctx = mm.Context(system, integrator, platform, {"Precision": "mixed"})
        ctx.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
        return ctx, integrator

    ctx, integrator = context(temp)
    ctx.setPositions(gro.getPositions())
    ctx.computeVirtualSites()
    mm.LocalEnergyMinimizer.minimize(ctx, 10., 2000)
    if phase == "liquid":  # melt with new velocities, then return to temp
        hot, hot_integrator = context(MELT_T)
        hot.setState(ctx.getState(getPositions=True))
        hot.setVelocitiesToTemperature(MELT_T*unit.kelvin, seed)
        hot_integrator.step(MELT_PS*500)
        ctx.setPositions(hot.getState(getPositions=True).getPositions())
        del hot
    ctx.setVelocitiesToTemperature(temp*unit.kelvin, seed)
    state = ctx.getState(getPositions=True)
    lattice = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)[::4].astype(float)
    box = np.diag(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom))
    equil_ps = EQUIL_PS["liquid" if phase == "liquid" else "crystal"]
    integrator.step(equil_ps*500)
    frames, energies, displacement = [], [], []
    for _ in range(FRAMES):
        integrator.step(STRIDE_PS*500)
        s = ctx.getState(getPositions=True, getEnergy=True)
        x = s.getPositions(asNumpy=True).value_in_unit(unit.angstrom).astype(np.float32)
        frames.append(x)
        energies.append(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        displacement.append(float(np.linalg.norm(minimum_image(x[::4]-lattice, box), axis=1).mean()))
    md_seconds = time.time() - start
    del ctx

    systemq, _ = system_ice.make_system(top)
    fq = fire.fire_integrator()
    cq = mm.Context(systemq, fq, platform, {"Precision": "mixed"})
    cq.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
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
    np.savez_compressed(OUT/f"{run}.npz", x=np.array(frames), xq=np.array(quenched), box=box,
                        energy=np.array(energies), displacement=np.array(displacement))
    simulated_ps = equil_ps + FRAMES*STRIDE_PS + (MELT_PS if phase == "liquid" else 0)
    info = dict(run=run, phase=phase, prep=prep, temperature=temp, velocity_seed=seed,
                structure=gro_path.name, equilibration_ps=equil_ps, frames=FRAMES, stride_ps=STRIDE_PS,
                simulated_ps=simulated_ps, md_seconds=md_seconds,
                ns_per_day=simulated_ps/1000/(md_seconds/86400), openmm=mm.__version__, **meta)
    (OUT/f"{run}.json").write_text(json.dumps(info, indent=1))
    print(f"{run}: {md_seconds:.0f} s, {info['ns_per_day']:.0f} ns/day, "
          f"energy {np.mean(energies[:10]):.0f} -> {np.mean(energies[-10:]):.0f} kJ/mol, "
          f"O displacement {displacement[0]:.2f} -> {displacement[-1]:.2f} A", flush=True)


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
    for run in (runs() if args.runs == ["all"] else args.runs):
        if (OUT/f"{run}.npz").exists() and not args.force:
            print(f"{run}: exists, skipped")
            continue
        simulate(run)


if __name__ == "__main__":
    main()
