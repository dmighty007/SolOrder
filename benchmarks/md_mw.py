"""Oxygen-only mW water (run D, model "mW"), in OpenMM without LAMMPS.

mW (Molinero & Moore 2009) is a Stillinger-Weber potential with
eps = 6.189 kcal/mol, sigma = 2.3925 A, a = 1.8, lambda = 23.15, gamma = 1.2,
cos(theta0) = -1/3, A = 7.049556277, B = 0.6022245584, p = 4, q = 0:
  E = sum_{i<j} A eps [B (sigma/r)^4 - 1] exp(sigma/(r - a sigma))
    + sum_i sum_{j<k} lambda eps (cos theta_jik + 1/3)^2
          exp(gamma sigma/(r_ij - a sigma)) exp(gamma sigma/(r_ik - a sigma))
The two-body term is a CustomNonbondedForce and the three-body term a
CustomManyParticleForce with a unique central particle; check_energy()
compares both with a direct numpy sum.

Systems: the run C oxygen positions of Ih, Ic and Isd (p1-p3) and the
seed-test liquid, at 220 K, which is also the TIP4P/Ice undercooling (mW
melts at 274.6 K). Protocol as md_conditions.py: NPT at 1 bar (anisotropic
for crystals; 200 ps, liquids melted at 300 K for 100 ps then 1 ns), NVT in
the final box for 50 ps, then 60 frames 5 ps apart, each FIRE-quenched for
300 steps. Time step 5 fs, Langevin 1/ps.

  python md_mw.py prepare          # locally: O-only inputs -> cache/conditions/mw_inputs/
  python md_mw.py check            # energy check against numpy
  python md_mw.py run <runs>       # MW_PLATFORM=CUDA|CPU, MW_DEVICE, MW_THREADS, MW_DIR
  python md_mw.py --list
Outputs: D_mW_<phase>_p<k>_220.npz/.json next to the inputs.
"""

import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
INPUTS = Path(os.environ.get("MW_DIR", HERE/"cache/conditions/mw_inputs"))
SEEDTEST = Path("/home/dm/Dibyendu/Projects/ice_nucleation_pg/data/seedtest")
sys.path[:0] = [str(HERE), str(SEEDTEST)]

KCAL = 4.184
EPS, SIGMA, A_CUT = 6.189*KCAL, .23925, 1.8  # kJ/mol, nm
LAMBDA, GAMMA, A_SW, B_SW = 23.15, 1.2, 7.049556277, .6022245584
RC = A_CUT*SIGMA
MASS = 18.015
DT_PS, FRICTION = .005, 1.
FRAMES, STRIDE_PS, NVT_PS, NPT_PS = 60, 5, 50, {"crystal": 200, "liquid": 1000}
MELT_PS, MELT_T, TEMP = 100, 300., 220.
PREPS = (1, 2, 3)


def runs():
    return [f"D_mW_{p}_p{k}_220" for p in ("Ih", "Ic", "Isd", "liquid") for k in PREPS]


def read_gro_oxygens(path):
    lines = Path(path).read_text().splitlines()
    atoms = lines[2:-1]
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in atoms if l[10:15].strip() == "OW"])
    box = [float(x) for x in lines[-1].split()]
    return xyz*10, np.diag(box[:3])*10


def write_gro(path, xyz, cell, title):
    lines = [title, str(len(xyz))]
    for m, p in enumerate(xyz/10):
        lines.append(f"{(m+1) % 100000:5d}{'SOL':<5s}{'OW':>5s}{(m+1) % 100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    c = cell/10
    lines.append(f"{c[0, 0]:.5f} {c[1, 1]:.5f} {c[2, 2]:.5f}")
    Path(path).write_text("\n".join(lines) + "\n")


def prepare():
    """O-only starting structures from run C (crystals) and the seed-test liquid."""
    replicas = HERE/"cache/replicas"
    INPUTS.mkdir(parents=True, exist_ok=True)
    for phase in ("Ih", "Ic", "Isd", "liquid"):
        for k in PREPS:
            source = SEEDTEST/"liquid.gro" if phase == "liquid" else replicas/f"{phase}_p{k}.gro"
            xyz, cell = read_gro_oxygens(source)
            write_gro(INPUTS/f"D_mW_{phase}_p{k}.gro", xyz, cell, f"mW O-only from {source.name}")
            stack = replicas/f"{phase}_p{k}.stack.npy"
            if stack.exists():
                np.save(INPUTS/f"D_mW_{phase}_p{k}.stack.npy", np.load(stack))
    for name in ("fire.py",):
        (INPUTS/name).write_text((SEEDTEST/name).read_text())
    (INPUTS/"md_mw.py").write_text(Path(__file__).read_text())
    print(f"inputs in {INPUTS}")


def make_system(n, cell_nm):
    import openmm as mm
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*[mm.Vec3(*v) for v in cell_nm])
    for _ in range(n):
        system.addParticle(MASS)
    two = mm.CustomNonbondedForce(
        "A*eps*(B*(sig/r)^4 - 1)*exp(sig/min(r - rc, -1e-6))")  # min: zero, not inf*0, past rc
    for name, value in (("A", A_SW), ("eps", EPS), ("B", B_SW), ("sig", SIGMA), ("rc", RC)):
        two.addGlobalParameter(name, value)
    two.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    two.setCutoffDistance(RC)
    two.setUseLongRangeCorrection(False)
    for _ in range(n):
        two.addParticle([])
    system.addForce(two)
    three = mm.CustomManyParticleForce(
        3, "lam*eps*(cos(angle(p2,p1,p3)) + 1/3)^2*exp(gam*sig/min(r12 - rc, -1e-6))"
           "*exp(gam*sig/min(r13 - rc, -1e-6)); r12 = distance(p1,p2); r13 = distance(p1,p3)")
    for name, value in (("lam", LAMBDA), ("eps", EPS), ("gam", GAMMA), ("sig", SIGMA), ("rc", RC)):
        three.addGlobalParameter(name, value)
    three.setPermutationMode(mm.CustomManyParticleForce.UniqueCentralParticle)
    three.setNonbondedMethod(mm.CustomManyParticleForce.CutoffPeriodic)
    three.setCutoffDistance(RC)
    for _ in range(n):
        three.addParticle([])
    system.addForce(three)
    return system


def numpy_energy(x, cell):
    """Direct mW energy (kJ/mol) for an orthorhombic cell, positions in nm."""
    box = np.diag(cell)
    d = x[None, :, :] - x[:, None, :]
    d -= box*np.round(d/box)
    r = np.linalg.norm(d, axis=2)
    np.fill_diagonal(r, np.inf)
    near = r < RC
    two = np.where(near, A_SW*EPS*(B_SW*(SIGMA/np.where(near, r, 1.))**4 - 1)
                   * np.exp(SIGMA/np.where(near, r - RC, -1.)), 0.)
    energy = two.sum()/2
    for i in range(len(x)):
        js = np.flatnonzero(near[i])
        u = d[i, js]/r[i, js, None]
        f = np.exp(GAMMA*SIGMA/(r[i, js] - RC))
        cos = u @ u.T
        jj, kk = np.triu_indices(len(js), 1)
        energy += np.sum(LAMBDA*EPS*(cos[jj, kk] + 1/3)**2*f[jj]*f[kk])
    return energy


def check_energy():
    import openmm as mm
    from openmm import unit
    xyz, cell = read_gro_oxygens(INPUTS/"D_mW_Ih_p1.gro")
    rng = np.random.default_rng(1)
    x = (xyz + rng.normal(0., .2, xyz.shape))/10  # nm; thermal-like disorder
    system = make_system(len(x), cell/10)
    ctx = mm.Context(system, mm.VerletIntegrator(.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(x)
    omm = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ref = numpy_energy(x, cell/10)
    per = omm/len(x)
    print(f"OpenMM {omm:.4f}  numpy {ref:.4f} kJ/mol  (relative difference {abs(omm-ref)/abs(ref):.2e}); "
          f"{per:.3f} kJ/mol per molecule")
    assert abs(omm - ref) < 1e-6*abs(ref) + 1e-3
    xi = xyz/10
    ctx.setPositions(xi)
    ideal = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)/len(xi)
    print(f"ideal Ih lattice at the run C spacing: {ideal:.3f} kJ/mol per molecule "
          f"({ideal/KCAL:.3f} kcal/mol; mW Ih minimum about -2 eps = {-2*EPS/KCAL:.3f})")


def platform_and_properties():
    import openmm as mm
    name = os.environ.get("MW_PLATFORM", "CPU")
    props = {}
    if name == "CPU" and "MW_THREADS" in os.environ:
        props["Threads"] = os.environ["MW_THREADS"]
    if name == "CUDA":
        props["Precision"] = "mixed"
        if "MW_DEVICE" in os.environ:
            props["DeviceIndex"] = os.environ["MW_DEVICE"]
    return mm.Platform.getPlatformByName(name), props


def simulate(run):
    import openmm as mm
    from openmm import unit
    import fire

    phase, prep = run.split("_")[2], int(run.split("_")[3][1:])
    xyz, cell = read_gro_oxygens(INPUTS/f"D_mW_{phase}_p{prep}.gro")
    n, liquid = len(xyz), phase == "liquid"
    seed = 1000*prep + sum(map(ord, run))
    platform, props = platform_and_properties()
    start = time.time()

    def context(t, barostat):
        system = make_system(n, cell/10)
        if barostat:
            if liquid:
                system.addForce(mm.MonteCarloBarostat(1.*unit.bar, t*unit.kelvin, 25))
            else:
                system.addForce(mm.MonteCarloAnisotropicBarostat(mm.Vec3(1., 1., 1.)*unit.bar, t*unit.kelvin,
                                                                 True, True, True, 25))
        integrator = mm.LangevinMiddleIntegrator(t*unit.kelvin, FRICTION/unit.picosecond, DT_PS*unit.picosecond)
        integrator.setRandomNumberSeed(seed)
        return mm.Context(system, integrator, platform, props), integrator

    def carry(src, dst):
        s = src.getState(getPositions=True)
        dst.setPeriodicBoxVectors(*s.getPeriodicBoxVectors())
        dst.setPositions(s.getPositions())

    steps_per_ps = round(1/DT_PS)

    def npt(ctx, integrator, ps):
        volumes = []
        for _ in range(ps):
            integrator.step(steps_per_ps)
            volumes.append(ctx.getState().getPeriodicBoxVolume().value_in_unit(unit.angstrom**3))
        return volumes

    ctx, integrator = context(TEMP, True)
    ctx.setPositions(xyz/10)
    mm.LocalEnergyMinimizer.minimize(ctx, 10., 2000)
    if liquid:
        hot, hot_integrator = context(MELT_T, True)
        carry(ctx, hot)
        hot.setVelocitiesToTemperature(MELT_T*unit.kelvin, seed)
        npt(hot, hot_integrator, MELT_PS)
        carry(hot, ctx)
        del hot
    ctx.setVelocitiesToTemperature(TEMP*unit.kelvin, seed)
    npt_ps = NPT_PS["liquid" if liquid else "crystal"]
    volumes = npt(ctx, integrator, npt_ps)
    state = ctx.getState(getPositions=True, getVelocities=True)
    del ctx

    nvt, integrator = context(TEMP, False)
    nvt.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    nvt.setPositions(state.getPositions())
    nvt.setVelocities(state.getVelocities())
    box = np.array(state.getPeriodicBoxVectors().value_in_unit(unit.angstrom))
    integrator.step(NVT_PS*steps_per_ps)
    frames, energies = [], []
    for _ in range(FRAMES):
        integrator.step(STRIDE_PS*steps_per_ps)
        s = nvt.getState(getPositions=True, getEnergy=True)
        frames.append(s.getPositions(asNumpy=True).value_in_unit(unit.angstrom).astype(np.float32))
        energies.append(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    md_seconds = time.time() - start
    del nvt
    lengths = np.diag(box)
    d = [x.astype(float) - frames[0] for x in frames]
    d = [v - lengths*np.round(v/lengths) for v in d]
    # No CMMotionRemover in this system: remove the whole-box drift before judging retention.
    displacement = [float(np.linalg.norm(v - v.mean(axis=0), axis=1).mean()) for v in d]

    systemq = make_system(n, box/10)
    fq = fire.fire_integrator()
    cq = mm.Context(systemq, fq, platform, props)
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
    np.savez_compressed(INPUTS/f"{run}.npz", x=np.array(frames), xq=np.array(quenched), cell=box,
                        energy=np.array(energies), displacement=np.array(displacement),
                        orientation=np.full(FRAMES, np.nan), npt_volume=np.array(volumes))
    simulated_ps = npt_ps + NVT_PS + FRAMES*STRIDE_PS + (MELT_PS if liquid else 0)
    out = dict(run=run, dataset="D", model="mW", phase=phase, prep=prep, temperature=TEMP, pressure_gpa=1e-4,
               structure=f"D_mW_{phase}_p{prep}.gro", velocity_seed=seed, molecules=n, npt_ps=npt_ps,
               nvt_ps=NVT_PS, frames=FRAMES, stride_ps=STRIDE_PS, timestep_ps=DT_PS, simulated_ps=simulated_ps,
               volume_start=volumes[0], volume_end=float(np.mean(volumes[-npt_ps//4:])),
               volume_drift_last_half=float(np.mean(volumes[-npt_ps//4:])/np.mean(volumes[npt_ps//2:3*npt_ps//4])-1),
               density_g_cm3=float(n*MASS/0.6022140857/abs(np.linalg.det(box))),
               md_seconds=md_seconds, ns_per_day=simulated_ps/1000/(md_seconds/86400),
               platform=platform.getName(), platform_properties=props, host=os.uname().nodename,
               openmm=mm.__version__)
    (INPUTS/f"{run}.json").write_text(json.dumps(out, indent=1))
    print(f"{run}: {md_seconds:.0f} s, {out['ns_per_day']:.1f} ns/day, density {out['density_g_cm3']:.3f}, "
          f"energy {np.mean(energies[:10])/n:.2f} -> {np.mean(energies[-10:])/n:.2f} kJ/mol/molecule, "
          f"O displacement {displacement[-1]:.2f} A", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", nargs="?", choices=("prepare", "check", "run"))
    parser.add_argument("runs", nargs="*")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    if args.list:
        print("\n".join(runs()))
    elif args.action == "prepare":
        prepare()
    elif args.action == "check":
        check_energy()
    elif args.action == "run":
        for run in (runs() if args.runs == ["all"] else args.runs):
            if (INPUTS/f"{run}.npz").exists():
                print(f"{run}: exists, skipped")
                continue
            simulate(run)


if __name__ == "__main__":
    main()
