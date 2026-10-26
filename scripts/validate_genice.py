"""Summarize SolOrder outputs on the checked-in, seeded GenIce fixtures.

Run regression assertions separately: python -m pytest -q tests/test_genice.py
"""

import argparse
import csv
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from SolOrder import OrderParameter, eFingerprint


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, default=root/"tests/data/genice")
    parser.add_argument("--output", type=Path, default=root/"docs/validation/genice_summary.csv")
    parser.add_argument("--n-jobs", type=int, default=2)
    args = parser.parse_args()
    manifest = json.loads((args.fixtures/"manifest.json").read_text())
    rows = []
    for phase, entry in manifest["structures"].items():
        u = mda.Universe(str(args.fixtures/entry["file"]))
        op = OrderParameter(u, cutoff=5., selection="name O", hydrogen_selection="name H", n_jobs=args.n_jobs)
        hb = op.hydrogen_bonds()
        values = {
            "cn_3p5": op.coordination_number(3.5),
            "smooth_cn_r0_3p5_cutoff_5": op.smooth_coordination(),
            "oto": op.OTO(), "tto": op.TTO(), "lsi_angstrom2": op.LSI(),
            "d5_angstrom": op.d5(), "f4_3p5": op.F4(cutoff=3.5),
            "hb_donated": hb.donated, "hb_accepted": hb.accepted,
            "q3_3p5": op.steinhardt(3, cutoff=3.5),
            "q4_3p5": op.steinhardt(4, cutoff=3.5),
            "q6_3p5": op.steinhardt(6, cutoff=3.5),
            "q3_averaged_3p5": op.steinhardt(3, averaged=True, cutoff=3.5),
            "q4_averaged_3p5": op.steinhardt(4, averaged=True, cutoff=3.5),
            "w6_normalized_3p5": op.steinhardt(6, invariant="w", normalized=True, cutoff=3.5),
            "connections_l6_threshold_0p5_cutoff_3p5": op.crystalline_connections(6, cutoff=3.5),
        }
        geometry = op.voronoi_geometry()
        values.update({
            "voronoi_volume_angstrom3": geometry.volume,
            "voronoi_area_angstrom2": geometry.area,
            "voronoi_density_angstrom_minus3": geometry.density,
            "voronoi_asphericity": geometry.asphericity,
            "centrosymmetry_n4_cutoff_5": op.centrosymmetry(4),
            "neighbor_distance_powers111_cutoff_5": op.neighbor_distance_parameter(),
            "angular_fourier_a1_unit_weight_cutoff_5": op.angular_fourier(),
            "bond_angle_n1_n1_phi0_cutoff_5": op.bond_angle_parameter(),
            "water_pair_angle_rad_cutoff_3p5": op.water_angular_distribution().angles,
        })
        neighborhood = op.neighborhood_parameters()
        histogram = op.angle_histogram_parameter()
        for channel in range(3):
            values[f"neighborhood_A{channel+1}_M1_cutoff_5"] = neighborhood[:, channel]
        for channel in range(4):
            values[f"angle_histogram_amplitude_f{channel+1}_bins180_cutoff_5"] = histogram[:, channel]
        classifications = {"chill": op.chill(), "chill_plus": op.chill_plus()}
        for local in (True, False):
            entropy = eFingerprint(u, selection="name O", local=local, n_jobs=args.n_jobs)
            name = "entropy_local" if local else "entropy_global"
            values[name] = entropy.Entropy()
            values[name+"_averaged"] = entropy.localEntropy()
        row = {"phase": phase, "waters": op.n, "density_angstrom_minus3": op.density()}
        for name, array in values.items():
            if not np.isfinite(array).all():
                raise ValueError(f"{phase}/{name} contains non-finite results")
            row[name+"_mean"] = float(np.mean(array))
            row[name+"_std"] = float(np.std(array))
        for name, result in classifications.items():
            for code, label in enumerate(("other", "Ih", "Ic", "ice_interface", "hydrate", "hydrate_interface")):
                row[f"{name}_{label}_count"] = int(np.count_nonzero(result.labels == code))
        rows.append(row)
        print(f"{phase:>4}: {op.n:4} waters; Q={row['oto_mean']:.6f}, "
              f"F4={row['f4_3p5_mean']:.6f}, HB={np.mean(hb.total):.1f}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
