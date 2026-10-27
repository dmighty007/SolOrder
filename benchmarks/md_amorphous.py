"""Amorphous ice pilot (TIP4P/Ice): HDA by pressure-induced amorphization, LDA two ways.

Legs, each continuing from the previous one where noted:
  compress_p<k>    ice Ih (run C structure) at 80 K, 0 -> 2.0 GPa
  decompress_p<k>  the compressed sample at 80 K, 2.0 GPa -> 1 bar (recovered HDA)
  heat_p<k>        recovered HDA at 1 bar, 80 -> 200 K (HDA -> LDA?)
  quench_p<k>      liquid at 300 K, 1 bar, cooled to 80 K (hyperquenched glass, LDA?)

Ramps are steps of 0.05 GPa or 5 K held for STEP_PS each (0.1 GPa/ns, 10 K/ns at
the default 500 ps), with an isotropic Monte Carlo barostat. Per step we record
the mean density over its second half; every SNAPSHOT_EVERY steps we save all
atom positions and the box. The pilot only checks densities and structure; it
does not yet produce benchmark data.

Needs OpenMM with CUDA (e.g. the gencompas conda environment):
  python md_amorphous.py compress_p1 decompress_p1 heat_p1 quench_p1
Outputs go to cache/amorphous/ (not committed): <leg>.npz, <leg>.json.
"""

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

import md_replicas

SEEDTEST = md_replicas.SEEDTEST
OUT = Path(__file__).resolve().parent/"cache/amorphous"
sys.path.insert(0, str(SEEDTEST))
STEP_PS, SNAPSHOT_EVERY, EQUIL_PS = 500, 5, 200
LOW_T, P_MAX, HEAT_T, LIQUID_T = 80., 2.0, 200., 300.
BAR_PER_GPA = 1e4


def ramp(start, stop, step):
    n = int(round(abs(stop-start)/step))
    return list(np.linspace(start, stop, n+1))


LEGS = {
    "compress": dict(source="Ih", control="pressure", values=ramp(0., P_MAX, .05), temperature=LOW_T),
    "decompress": dict(source="compress", control="pressure", values=ramp(P_MAX, 0., .05), temperature=LOW_T),
    "heat": dict(source="decompress", control="temperature", values=ramp(LOW_T, HEAT_T, 5.), pressure=0.),
    "quench": dict(source="liquid", control="temperature", values=ramp(LIQUID_T, LOW_T, 5.), pressure=0.),
}


def simulate(name):
    import openmm as mm
    from openmm import app, unit
    import system_ice

    leg, prep = name.split("_")
    spec = LEGS[leg]
    prep_index = int(prep[1:])
    seed = 4000 + 100*prep_index + list(LEGS).index(leg)
    platform = mm.Platform.getPlatformByName("CUDA")
    source = spec["source"]
    fresh = source in ("Ih", "liquid")
    gro_path, _ = md_replicas.prepare("liquid" if leg == "quench" else "Ih", prep_index)
    gro = app.GromacsGroFile(str(gro_path))
    top = app.GromacsTopFile(str(gro_path.with_suffix(".top")), periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir=str(SEEDTEST))
    system, integrator = system_ice.make_system(top, temp=spec.get("temperature", spec["values"][0]))
    barostat = mm.MonteCarloBarostat(1.*unit.bar, LOW_T*unit.kelvin, 25)
    system.addForce(barostat)
    integrator.setRandomNumberSeed(seed)
    ctx = mm.Context(system, integrator, platform, {"Precision": "mixed"})
    n = system.getNumParticles()//4
    mass = n*18.01528/0.6022140857  # g/cm3 times A^3

    def set_conditions(temperature, pressure):
        integrator.setTemperature(temperature*unit.kelvin)
        ctx.setParameter(mm.MonteCarloBarostat.Pressure(), pressure*BAR_PER_GPA)
        ctx.setParameter(mm.MonteCarloBarostat.Temperature(), temperature)

    start = time.time()
    if fresh:  # minimize and equilibrate at the first condition (the liquid starts at 300 K)
        ctx.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
        ctx.setPositions(gro.getPositions())
        ctx.computeVirtualSites()
        mm.LocalEnergyMinimizer.minimize(ctx, 10., 2000)
        t0 = spec.get("temperature", spec["values"][0])
        p0 = spec.get("pressure", spec["values"][0])
        set_conditions(t0, p0)
        ctx.setVelocitiesToTemperature(t0*unit.kelvin, seed)
        integrator.step(EQUIL_PS*500)
    else:  # continue from the end of the source leg
        previous = np.load(OUT/f"{source}_{prep}.npz")
        ctx.setPeriodicBoxVectors(*(previous["final_cell"]/10.))
        ctx.setPositions(previous["final_x"]/10.)
        ctx.setVelocities(previous["final_v"])
    density, snapshots, cells, controls = [], [], [], []
    for k, value in enumerate(spec["values"]):
        temperature = value if spec["control"] == "temperature" else spec["temperature"]
        pressure = value if spec["control"] == "pressure" else spec["pressure"]
        set_conditions(temperature, pressure)
        volumes = []
        for _ in range(STEP_PS):
            integrator.step(500)
            volumes.append(ctx.getState().getPeriodicBoxVolume().value_in_unit(unit.angstrom**3))
        density.append(mass/np.mean(volumes[STEP_PS//2:]))
        if k % SNAPSHOT_EVERY == 0 or k == len(spec["values"])-1:
            s = ctx.getState(getPositions=True)
            snapshots.append(s.getPositions(asNumpy=True).value_in_unit(unit.angstrom).astype(np.float32))
            cells.append(np.array(s.getPeriodicBoxVectors().value_in_unit(unit.angstrom)))
            controls.append(value)
        print(f"{name} {spec['control']} {value:.2f}: density {density[-1]:.4f} g/cm3, "
              f"{time.time()-start:.0f} s", flush=True)
    final = ctx.getState(getPositions=True, getVelocities=True)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT/f"{name}.npz", control=np.array(spec["values"]), density=np.array(density),
                        snapshot_control=np.array(controls), x=np.array(snapshots), cell=np.array(cells),
                        final_x=final.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
                        final_v=final.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond),
                        final_cell=np.array(final.getPeriodicBoxVectors().value_in_unit(unit.angstrom)))
    info = dict(leg=name, **{k: v for k, v in spec.items() if k != "values"}, start=spec["values"][0],
                stop=spec["values"][-1], step_ps=STEP_PS, steps=len(spec["values"]), molecules=n,
                structure=gro_path.name, velocity_seed=seed, md_seconds=time.time()-start,
                openmm=mm.__version__)
    (OUT/f"{name}.json").write_text(json.dumps(info, indent=1))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("legs", nargs="+")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for name in args.legs:
        simulate(name)


if __name__ == "__main__":
    main()
