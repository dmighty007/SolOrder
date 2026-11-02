"""Convergence of the fixed 300-step FIRE quench used for every "quenched" frame.

For representative thermal frames (TIP4P/Ice Ih, Ic and liquid from run C1;
ices V and VIII and liquid at 1.0 GPa from run E; TIP4P/2005 Ih and liquid
from run D; mW Ih and liquid), FIRE is continued from the thermal frame with
the benchmark settings and checkpointed after 100, 300, 1000, 3000 and 10000
steps. Step "quench" records the potential energy and RMS net force per
molecule (kJ/mol/nm) at each checkpoint and saves the positions (OpenMM environment); step "compare"
computes the registry parameters at 300 and 10000 steps and reports their
change relative to the Ih - liquid separation of each parameter at 220 K
(SolOrder environment).

  <gencompas>/bin/python quench_check.py quench    # cache/quench_check.npz
  ../.conda-env/bin/python quench_check.py compare  # results/quench_check.csv, quench_convergence.csv
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
CHECKPOINTS = (0, 100, 300, 1000, 3000, 10000)
SAMPLES = {  # name: (npz, frame, model, gro for topology)
    "TIP4P/Ice Ih": ("replicas/Ih_p1_220.npz", 0, "tip4pice", "replicas/Ih_p1.gro"),
    "TIP4P/Ice Ic": ("replicas/Ic_p1_220.npz", 0, "tip4pice", "replicas/Ic_p1.gro"),
    "TIP4P/Ice liquid": ("replicas/liquid_p1_220.npz", 0, "tip4pice", "replicas/liquid_p1.gro"),
    "TIP4P/Ice V 0.5 GPa": ("conditions/E_V_p1_220_0.5GPa.npz", 0, "tip4pice", "conditions/E_V_p1.gro"),
    "TIP4P/Ice VIII 3 GPa": ("conditions/E_VIII_p1_220_3.0GPa.npz", 0, "tip4pice", "conditions/E_VIII_p1.gro"),
    "TIP4P/Ice liquid 1 GPa": ("conditions/E_liquid_p1_220_1.0GPa.npz", 0, "tip4pice", "replicas/liquid_p1.gro"),
    "TIP4P/2005 Ih": ("conditions/D_tip4p2005_Ih_p1_220.npz", 0, "tip4p2005", "replicas/Ih_p1.gro"),
    "TIP4P/2005 liquid": ("conditions/D_tip4p2005_liquid_p1_220.npz", 0, "tip4p2005", "replicas/liquid_p1.gro"),
    "mW Ih": ("conditions/mw_inputs/D_mW_Ih_p1_220.npz", 0, "mW", "conditions/mw_inputs/D_mW_Ih_p1.gro"),
    "mW liquid": ("conditions/mw_inputs/D_mW_liquid_p1_220.npz", 0, "mW",
                  "conditions/mw_inputs/D_mW_liquid_p1.gro"),
}


def quench():
    import openmm as mm
    from openmm import app, unit
    import md_conditions
    import md_mw
    sys.path.insert(0, str(md_conditions.SEEDTEST))
    import fire
    import system_ice

    platform = mm.Platform.getPlatformByName("CUDA")
    saved, rows = {}, []
    for name, (npz, frame, model, gro_path) in SAMPLES.items():
        data = np.load(CACHE/npz)
        cell = data["cell"] if "cell" in data.files else np.diag(data["box"])
        x = data["x"][frame].astype(float)
        if model == "mW":
            system = md_mw.make_system(len(x), cell/10)
        else:
            gro = app.GromacsGroFile(str(CACHE/gro_path))
            top = app.GromacsTopFile(str((CACHE/gro_path).with_suffix(".top")),
                                     periodicBoxVectors=gro.getPeriodicBoxVectors(),
                                     includeDir=str(md_conditions.SEEDTEST))
            system, _ = system_ice.make_system(top)
            md_conditions.set_model(system, model)
        integrator = fire.fire_integrator()
        ctx = mm.Context(system, integrator, platform, {"Precision": "mixed"})
        ctx.setPeriodicBoxVectors(*(cell/10.))
        fire.set_inverse_masses(integrator, system)
        ctx.setPositions(x/10.)
        ctx.setVelocitiesToTemperature(0.)
        integrator.setGlobalVariableByName("dtc", .001)
        integrator.setGlobalVariableByName("alpha", .1)
        integrator.setGlobalVariableByName("npos", 0)
        done = 0
        n = len(x)//(1 if model == "mW" else 4)
        for checkpoint in CHECKPOINTS:
            integrator.step(checkpoint - done)
            done = checkpoint
            state = ctx.getState(getPositions=True, getEnergy=True, getForces=True)
            forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
            # Net force per molecule over O, H, H: constraint forces of rigid water cancel within a
            # molecule, and OpenMM reports the M-site force in addition to its redistribution.
            net = forces.reshape(n, -1, 3)[:, :3].sum(axis=1)
            rows.append(dict(sample=name, steps=checkpoint,
                             energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)/n,
                             rms_force=float(np.sqrt(np.mean(np.sum(net**2, axis=1))))))
            saved[f"{name}|{checkpoint}"] = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        saved[f"{name}|cell"] = cell
        print(name, " ".join(f"{r['steps']}: {r['energy']:.3f} ({r['rms_force']:.0f})" for r in rows[-len(CHECKPOINTS):]),
              flush=True)
    np.savez_compressed(CACHE/"quench_check.npz", **saved)
    (CACHE/"quench_check.json").write_text(json.dumps(rows, indent=1))


def compare():
    import MDAnalysis as mda
    import pandas as pd
    from compute import evaluate, triclinic_box

    rows = json.loads((CACHE/"quench_check.json").read_text())
    convergence = pd.DataFrame(rows)
    final = convergence.groupby("sample").energy.transform("last")
    convergence["energy_above_final"] = convergence.energy - final
    RESULTS.mkdir(exist_ok=True)
    convergence.to_csv(RESULTS/"quench_convergence.csv", index=False)
    data = np.load(CACHE/"quench_check.npz")
    prep = pd.read_pickle(CACHE/"prep.pkl")
    ref = prep[(prep.kind == "thermal") & (prep.temperature == 220.)]
    out, all_values = [], {}
    for name, (_, _, model, gro_path) in SAMPLES.items():
        values = all_values[name] = {}
        for steps in (300, CHECKPOINTS[-1]):
            x = data[f"{name}|{steps}"]
            u = mda.Universe(str(CACHE/gro_path))
            hydrogen = None if model == "mW" else "name HW1 HW2"
            u.atoms.positions = x
            u.dimensions = triclinic_box(*data[f"{name}|cell"])
            values[steps] = evaluate(u, "name OW", hydrogen)
        for column, a in values[300].items():
            b = values[CHECKPOINTS[-1]][column]
            if column not in ref or not np.isfinite(a).any():
                continue
            scale = abs(ref.loc[ref.phase == "Ih", column].median() - ref.loc[ref.phase == "liquid", column].median())
            finite = np.isfinite(a) & np.isfinite(b)
            out.append(dict(sample=name, column=column, mean_shift=float(np.mean(b[finite]-a[finite])),
                            mean_abs_change=float(np.mean(np.abs(b[finite]-a[finite]))), ih_liquid_scale=scale,
                            relative=float(np.mean(np.abs(b[finite]-a[finite]))/scale) if scale > 0 else np.nan))
    table = pd.DataFrame(out)
    table.to_csv(RESULTS/"quench_check.csv", index=False)
    # Effect on the scores: Ih versus liquid oracle score with both at 300 steps or both at 10000.
    from metrics import ks
    scores = []
    for ice, liquid in (("TIP4P/Ice Ih", "TIP4P/Ice liquid"), ("TIP4P/Ice Ih", "TIP4P/Ice Ic"),
                        ("TIP4P/2005 Ih", "TIP4P/2005 liquid"), ("mW Ih", "mW liquid")):
        for column in all_values[ice][300]:
            a300, b300 = all_values[ice][300][column], all_values[liquid][300][column]
            a_end, b_end = all_values[ice][CHECKPOINTS[-1]][column], all_values[liquid][CHECKPOINTS[-1]][column]
            if not np.isfinite(a300).any() or column in ("chill", "chill_plus", "ptm"):
                continue
            scores.append(dict(pair=f"{ice} / {liquid}", column=column, oracle_300=(1+ks(a300, b300))/2,
                               oracle_converged=(1+ks(a_end, b_end))/2))
    scores = pd.DataFrame(scores)
    scores["change"] = scores.oracle_converged - scores.oracle_300
    scores.to_csv(RESULTS/"quench_scores.csv", index=False)
    print(scores.groupby("pair").change.agg(["median", lambda c: c.abs().max()]).round(3).to_string())
    print(scores.loc[scores.change.abs() > .02].round(3).to_string(index=False))
    pd.set_option("display.width", 220)
    print(convergence.pivot(index="sample", columns="steps", values="energy_above_final").round(4).to_string())
    print(convergence.pivot(index="sample", columns="steps", values="rms_force").round(0).to_string())
    print(table.pivot(index="column", columns="sample", values="relative").round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("quench", "compare"))
    quench() if parser.parse_args().step == "quench" else compare()
