"""Independent numerical spot checks on representative nucleation snapshots."""

import argparse
import json
from pathlib import Path

import freud
import numpy as np
from MDAnalysis.lib.distances import calc_dihedrals, calc_angles, minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors
from scipy.special import eval_legendre

from SolOrder import OrderParameter
from SolOrder.trajectory_io import snapshot_template
from validate_nucleation import iter_frames, oxygen_universe, write_csv


AUDIT_SOURCES = {"PG_E", "PG_I_Ic", "PG_M_Ih", "TIP4P_w0", "TIP4P_w1", "TIP4P_w2",
                 "mW_brute_T200", "mW_brute_T205", "mW_brute_T210",
                 "PG_I_quenched_movie", "PG_M_quenched_movie"}


def compare(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        raise AssertionError("Reference NaN pattern differs")
    finite = np.isfinite(a) & np.isfinite(b)
    return float(np.max(np.abs(a[finite]-b[finite]))) if finite.any() else 0.


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=root/"docs/validation/nucleation")
    args = parser.parse_args()
    manifest = json.loads((args.directory/"manifest.json").read_text())
    rows = []
    for source in manifest["sources"]:
        if source["name"] not in AUDIT_SOURCES:
            continue
        source = dict(source, selected_frames=sorted(set([source["selected_frames"][0], source["selected_frames"][-1]])))
        if "boxes" in source:
            source["boxes"] = {int(k): v for k, v in source["boxes"].items()}
        u = snapshot_template(source["topology"]) if source["full_atom"] else oxygen_universe(source["atoms"])
        for index, positions, dimensions, _ in iter_frames(source):
            u.atoms.positions = positions
            u.dimensions = dimensions
            op = OrderParameter(u, cutoff=5, selection="name O OW", hydrogen_selection="name H H1 H2 HW1 HW2")
            np.testing.assert_array_equal(op.pos, np.asarray(positions, dtype=np.float32)[op.atom_indices])
            row = {"source": source["name"], "frame": index}
            graph = op.neighbor_list(3.5)
            source_rows = np.repeat(np.arange(op.n), graph.counts)
            nlist = freud.locality.NeighborList.from_arrays(op.n, op.n, source_rows,
                                                           graph.indices, graph.vectors)
            box = freud.box.Box.from_matrix(triclinic_vectors(op.box, dtype=float).T)
            center_rows = np.linspace(0, op.n-1, 32, dtype=int)
            for l in (3, 4, 6):
                direct = []
                for i in center_rows:
                    sl = slice(graph.offsets[i], graph.offsets[i+1])
                    if graph.counts[i] == 0:
                        direct.append(np.nan)
                        continue
                    unit = graph.vectors[sl]/graph.distances[sl, None]
                    q_squared = eval_legendre(l, np.clip(unit @ unit.T, -1, 1)).mean()
                    direct.append(np.sqrt(max(0., q_squared)))
                row[f"q{l}_Legendre_max_error"] = compare(op.steinhardt(l, cutoff=3.5)[center_rows], direct)
                for averaged in (False, True):
                    ref = freud.order.Steinhardt(l=l, average=averaged).compute((box, op.pos), neighbors=nlist)
                    name = f"q{l}"+("_averaged" if averaged else "")
                    row[name+"_max_error"] = compare(op.steinhardt(l, averaged=averaged, cutoff=3.5), ref.particle_order)
            for normalized in (False, True):
                ref = freud.order.Steinhardt(l=6, wl=True, wl_normalize=normalized).compute((box, op.pos), neighbors=nlist)
                row[("w6_normalized" if normalized else "w6")+"_max_error"] = compare(
                    op.steinhardt(6, invariant="w", normalized=normalized, cutoff=3.5), ref.particle_order)
            g = op.neighbor_list(5.)
            centers = np.linspace(0, op.n-1, 32, dtype=int)
            q_ref, tto_ref, lsi_ref = [], [], []
            ambiguous = []
            for i in centers:
                radii = g.distances[g.offsets[i]:g.offsets[i+1]]
                vectors = g.vectors[g.offsets[i]:g.offsets[i+1]]
                if len(radii) >= 4:
                    unit = vectors[:4]/radii[:4, None]
                    cosine = unit @ unit.T
                    q_ref.append(1-3/8*np.sum((cosine[np.triu_indices(4, 1)]+1/3)**2))
                    tto_ref.append(1-np.sum((radii[:4]-radii[:4].mean())**2)/(12*radii[:4].mean()**2))
                else:
                    q_ref.append(np.nan); tto_ref.append(np.nan)
                n = int(np.sum(radii < 3.7))
                lsi_ref.append(float(np.var(np.diff(radii[:n+1]))) if n > 0 and n < len(radii) else np.nan)
            row["Q_max_error"] = compare(op.OTO()[centers], q_ref)
            row["TTO_max_error"] = compare(op.TTO()[centers], tto_ref)
            row["LSI_max_error"] = compare(op.LSI()[centers], lsi_ref)
            if source["full_atom"]:
                oh = op._water_geometry()
                oo = graph.vectors
                ha, hb = oh[source_rows], oh[graph.indices]
                ia = np.argmax(np.linalg.norm(ha-oo[:, None, :], axis=2), axis=1)
                ib = np.argmax(np.linalg.norm(hb+oo[:, None, :], axis=2), axis=1)
                a = ha[np.arange(len(oo)), ia]
                b = oo+hb[np.arange(len(oo)), ib]
                angle = calc_dihedrals(a, np.zeros_like(oo), oo, b)
                ref = np.bincount(source_rows, weights=np.cos(3*angle), minlength=op.n)/graph.counts
                row["F4_max_error"] = compare(op.F4(cutoff=3.5), ref)
                distribution = op.water_angular_distribution()
                selected = np.linspace(0, len(distribution.pairs)-1, 64, dtype=int)
                i, j = distribution.pairs[selected].T
                delta = minimize_vectors(op.pos[j]-op.pos[i], op.box)
                zero = np.zeros_like(delta)
                angles = [calc_angles(oh[i, h], zero, delta) for h in range(2)]
                angles += [calc_angles(delta+oh[j, h], delta, zero) for h in range(2)]
                row["Ptheta_max_error"] = compare(distribution.angles[selected], np.min(angles, axis=0))
            errors = [v for k, v in row.items() if k.endswith("_max_error")]
            passed = []
            for key, value in row.items():
                if not key.endswith("_max_error"):
                    continue
                tolerance = 1e-7 if "Legendre" in key else 1e-10 if key in ("Q_max_error", "TTO_max_error", "LSI_max_error") else 1e-6 if key in ("F4_max_error", "Ptheta_max_error") else 1e-4
                passed.append(value <= tolerance)
            row["passed"] = all(passed)
            rows.append(row)
            print(source["name"],index,"max_error",max(errors),"passed",row["passed"],flush=True)
            write_csv(args.directory/"independent_audits.csv", rows)
    if not all(r["passed"] for r in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
