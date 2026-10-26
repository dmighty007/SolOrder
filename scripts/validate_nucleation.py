"""Validate SolOrder on local TIP4P/Ice, PathGennie and mW trajectories.

Inputs are read-only. Results and selected-snapshot hashes are written under
SolOrder. Large source trajectories are never copied into the repository.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDReader
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import rankdata
from importlib.metadata import version

from SolOrder import OrderParameter, eFingerprint
from SolOrder.trajectory_io import sample_indices, scan_coordinate_chunks, read_coordinate_chunk, index_xyz, read_xyz_frame, snapshot_template


NVT_BOX = 50.5301  # A; inputs/start.gro and the project's movie/analysis code.
LABELS = ("other", "Ih", "Ic", "ice_interface", "hydrate", "hydrate_interface")


def oxygen_universe(n):
    u = mda.Universe.empty(n, n_residues=n, atom_resindex=np.arange(n), trajectory=True)
    u.add_TopologyAttr("names", ["OW"]*n)
    u.add_TopologyAttr("resids", np.arange(1, n+1))
    return u


def write_csv(path, rows):
    if not rows:
        return
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def largest_core(op, labels, include=(1, 2)):
    selected = np.isin(labels, include)
    if not selected.any():
        return np.zeros(op.n, dtype=bool)
    g = op.neighbor_list(3.5)
    source = np.repeat(np.arange(op.n), g.counts)
    mask = selected[source] & selected[g.indices]
    adjacency = coo_matrix((np.ones(mask.sum()), (source[mask], g.indices[mask])), shape=(op.n, op.n))
    _, components = connected_components(adjacency, directed=False)
    counts = np.bincount(components[selected])
    return selected & (components == counts.argmax())


def fields(op, full_atom):
    values = {
        "cn_3p5": op.coordination_number(3.5), "smooth_cn": op.smooth_coordination(),
        "Q": op.OTO(), "Q_averaged": op.tetrahedral_order(averaged=True),
        "TTO": op.TTO(), "LSI": op.LSI(), "d5": op.d5(),
        "q3": op.steinhardt(3, cutoff=3.5), "q4": op.steinhardt(4, cutoff=3.5),
        "q6": op.steinhardt(6, cutoff=3.5),
        "q3_averaged": op.steinhardt(3, averaged=True, cutoff=3.5),
        "q4_averaged": op.steinhardt(4, averaged=True, cutoff=3.5),
        "q6_averaged": op.steinhardt(6, averaged=True, cutoff=3.5),
        "w6": op.steinhardt(6, invariant="w", cutoff=3.5),
        "w6_normalized": op.steinhardt(6, invariant="w", normalized=True, cutoff=3.5),
        "connections6": op.crystalline_connections(6, cutoff=3.5),
        "CSP4": op.centrosymmetry(4), "D": op.neighbor_distance_parameter(),
        "F": op.angular_fourier(), "B": op.bond_angle_parameter(),
    }
    voronoi = op.voronoi_geometry()
    values.update(voronoi_volume=voronoi.volume, voronoi_area=voronoi.area,
                  voronoi_density=voronoi.density, asphericity=voronoi.asphericity)
    neighborhood, histogram = op.neighborhood_parameters(), op.angle_histogram_parameter()
    for i in range(3):
        values[f"A{i+1}"] = neighborhood[:, i]
    for i in range(4):
        values[f"H{i+1}"] = histogram[:, i]
    for local, name in ((True, "entropy_local"), (False, "entropy_global")):
        entropy = eFingerprint(op.u, selection=op.selection, local=local, n_jobs=op.n_jobs, backend=op.backend)
        values[name] = entropy.Entropy()
        values[name+"_averaged"] = entropy.localEntropy()
    pair_angles = None
    if full_atom:
        hb = op.hydrogen_bonds()
        values.update(F4=op.F4(cutoff=3.5), HB=hb.total, HB_donated=hb.donated, HB_accepted=hb.accepted)
        pair_angles = op.water_angular_distribution()
        # TIP4P/Ice inputs define H=+0.5897e, M=-1.1794e. Derive the
        # actual dipole and M-site offset from this frame, not a fitted value.
        m_sites = op.u.select_atoms("name M MW")
        by_residue = {a.resindex: a.position.astype(float) for a in m_sites}
        if len(by_residue) == op.n:
            om = np.array([by_residue[r] for r in op._resindices])-op.pos
            from MDAnalysis.lib.distances import minimize_vectors
            om = minimize_vectors(om, op.dimensions)
            oh = op._water_geometry()
            mu = np.linalg.norm(.5897*(oh[:, 0]+oh[:, 1])-1.1794*om, axis=1)*4.80320471257
            values["Q_T_Debye_A"] = op.quadrupole_parameter(mu, negative_site_distance=np.linalg.norm(om, axis=1))
    return values, pair_angles


def summarize(array):
    array = np.asarray(array)
    finite = array[np.isfinite(array)]
    return {
        "mean": float(finite.mean()) if len(finite) else None,
        "std": float(finite.std()) if len(finite) else None,
        "nan_count": int(np.isnan(array).sum()), "inf_count": int(np.isinf(array).sum()),
        "size": int(array.size),
    }


def core_auc(array, core):
    finite = np.isfinite(array)
    values, truth = np.asarray(array)[finite], core[finite]
    positive, negative = int(truth.sum()), int((~truth).sum())
    if not positive or not negative:
        return None
    ranks = rankdata(values)
    return float((ranks[truth].sum()-positive*(positive+1)/2)/(positive*negative))


def bound_violations(values):
    bounds = {"Q": (-3, 1), "Q_averaged": (-3, 1), "TTO": (0, 1),
              "q3": (0, 1), "q4": (0, 1), "q6": (0, 1),
              "q3_averaged": (0, 1), "q4_averaged": (0, 1), "q6_averaged": (0, 1),
              "F4": (-1, 1), "F": (-1, 1), "B": (-1, 1),
              "LSI": (0, np.inf), "CSP4": (0, np.inf),
              "asphericity": (1, np.inf), "voronoi_volume": (0, np.inf)}
    for i in range(1, 5):
        bounds[f"H{i}"] = (0, 1)
    for name in ("entropy_local", "entropy_global", "entropy_local_averaged", "entropy_global_averaged"):
        bounds[name] = (-np.inf, 0)
    count = 0
    for name, (lower, upper) in bounds.items():
        if name in values:
            a = values[name]
            count += int(np.count_nonzero((a < lower-1e-8) | (a > upper+1e-8)))
    return count


def reference_hb(op, reference):
    # Separate NumPy/cKDTree implementation of the same donor-H-acceptor count.
    pairs, delta = reference.neighbor_pairs(op.pos/10, op.dimensions[:3]/10)
    source = np.r_[pairs[:, 0], pairs[:, 1]]
    target = np.r_[pairs[:, 1], pairs[:, 0]]
    vectors = np.concatenate((delta, -delta))
    oh = op._water_geometry()[source]/10
    cosine = np.sum(oh*vectors[:, None, :], axis=2)/(np.linalg.norm(oh, axis=2)*np.linalg.norm(vectors, axis=1)[:, None])
    donor_edge, _ = np.nonzero(cosine > np.cos(np.pi/6))
    return (np.bincount(source[donor_edge], minlength=op.n)+np.bincount(target[donor_edge], minlength=op.n))


def prepare(root, output, frames, shooting_frames, include_quenched):
    """Return small source manifests and lazy frame iterators; no data mutation."""
    import h5py
    sources = []
    for tag, name in (("E", "ck_E.h5"), ("F", "ck_F.h5"), ("G_Ic", "ck_G.h5"),
                      ("I_Ic", "ck_I_final.h5"), ("J_Ih", "ck_J.h5"), ("K_Ih", "ck_K.h5"),
                      ("M_Ih", "ck_M_final.h5"), ("O_Ic", "ck_O_Ic_final.h5")):
        path = root/"data/pg_checkpoints"/name
        if not path.exists():
            continue
        entry = dict(name="PG_"+tag, family="PathGennie_raw", path=str(path), kind="hdf5",
                     topology=str(root/"inputs/start.gro"), full_atom=True,
                     box_provenance="Fixed NVT box from inputs/start.gro and project analysis code", box_approximate=False)
        try:
            with h5py.File(path, "r") as handle:
                entry["frames"] = len(handle["trajectory"])
                entry["atoms"] = handle["trajectory"].shape[1]
        except OSError as exc:
            entry.update(kind="zlib_recovery", header_error=str(exc), atoms=16384)
            print(f"Recovering coordinate index: {name}", flush=True)
            records = scan_coordinate_chunks(path, entry["atoms"])
            entry["frames"] = len(records)
            entry["recovery_records"] = records
            entry["frame_order"] = "Physical zlib chunk order; logical HDF5 ordering unavailable"
            movie_tag = tag.split("_")[0]
            info_path = root/f"outputs/pg_{movie_tag}_movie/pg_movie_info_all.dat"
            if info_path.exists():
                entry["source_movie_reported_frames"] = len(np.loadtxt(info_path))
                entry["unrecovered_relative_to_movie"] = entry["source_movie_reported_frames"]-len(records)
            print(f"Recovered {len(records)} coordinate chunks: {name}", flush=True)
        if tag == "O_Ic":
            info_path = root/"outputs/pg_O_movie/pg_movie_info_all.dat"
            info = np.loadtxt(info_path)
            boxes = {int(row[0]): float(row[-1]) for row in info}
            if len(boxes) != entry["frames"] or set(boxes) != set(range(entry["frames"])):
                raise ValueError("NPT reconstructed box rows do not cover the checkpoint")
            entry.update(boxes=boxes, box_approximate=True,
                         box_provenance=str(info_path)+"; source project's cycle-window median reconstruction")
        entry["selected_frames"] = sample_indices(entry["frames"], frames)
        sources.append(entry)
    for path in sorted((root/"outputs/tip4pice_T220").rglob("traj.dcd")):
        relative = path.parent.relative_to(root/"outputs/tip4pice_T220")
        with DCDReader(str(path)) as reader:
            count, atoms = len(reader), reader.n_atoms
        entry = dict(name="TIP4P_"+str(relative).replace("/", "_"), family="MD_shooting" if len(relative.parts)>1 else "metadynamics",
                     path=str(path), kind="dcd", topology=str(root/"inputs/w1000.pdb"),
                     frames=count, atoms=atoms, full_atom=True, box_provenance="Native per-frame DCD unit cell", box_approximate=False)
        entry["selected_frames"] = sample_indices(count, shooting_frames if entry["family"] == "MD_shooting" else frames)
        sources.append(entry)
    for path in sorted((root/"outputs/mw_nucleation").rglob("traj.npz")):
        with np.load(path) as archive:
            count = len(archive["t"])
            atoms = archive["pos"].shape[1]
        entry = dict(name="mW_"+path.parent.name, family="mW_bruteforce", path=str(path), kind="npz",
                     frames=count, atoms=atoms, full_atom=False, box_provenance="Native per-frame NPZ box (nm converted to A)", box_approximate=False,
                     selected_frames=sample_indices(count, frames))
        sources.append(entry)
    if include_quenched:
        for tag in ("E", "I", "M", "O"):
            path = root/f"outputs/pg_{tag}_movie/pg_movie.xyz"
            records = index_xyz(path)
            entry = dict(name="PG_"+tag+"_quenched_movie", family="PathGennie_quenched", path=str(path), kind="xyz",
                         frames=len(records), atoms=records[0][1], full_atom=False,
                         box_provenance="XYZ comment (reconstructed NPT box)" if tag == "O" else "Fixed NVT box from source project",
                         box_approximate=tag == "O", xyz_records=records, selected_frames=sample_indices(len(records), frames))
            sources.append(entry)
    for entry in sources:
        stat = Path(entry["path"]).stat()
        entry.update(size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    return sources


def iter_frames(source):
    import h5py
    path = source["path"]
    if source["kind"] == "hdf5":
        with h5py.File(path, "r") as archive:
            for i in source["selected_frames"]:
                length = source.get("boxes", {}).get(i, NVT_BOX)
                yield i, archive["trajectory"][i], np.r_[[length]*3, [90.]*3], None
    elif source["kind"] == "zlib_recovery":
        for i in source["selected_frames"]:
            yield i, read_coordinate_chunk(path, source["recovery_records"][i], source["atoms"]), np.r_[[NVT_BOX]*3, [90.]*3], None
    elif source["kind"] == "dcd":
        with DCDReader(path) as reader:
            for i in source["selected_frames"]:
                ts = reader[i]
                if ts.dimensions is None:
                    raise ValueError("DCD frame has no native periodic box")
                yield i, ts.positions.copy(), ts.dimensions.copy(), float(ts.time)/1000
    elif source["kind"] == "npz":
        with np.load(path) as archive:
            positions, boxes, times = archive["pos"], archive["box"], archive["t"]
            for i in source["selected_frames"]:
                yield i, positions[i]*10., np.r_[boxes[i]*10., [90.]*3], float(times[i])
    else:
        for i in source["selected_frames"]:
            positions, length, original = read_xyz_frame(path, source["xyz_records"][i])
            yield i, positions, np.r_[[NVT_BOX if length is None else length]*3, [90.]*3], None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--data-root", type=Path, default=Path("/home/dm/Dibyendu/Projects/ice_nucleation_pg"))
    parser.add_argument("--output", type=Path, default=root/"docs/validation/nucleation")
    parser.add_argument("--frames", type=int, default=12)
    parser.add_argument("--shooting-frames", type=int, default=3)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--no-quenched", action="store_true")
    parser.add_argument("--only", help="Comma-separated source names")
    parser.add_argument("--limit-sources", type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("nucleation_reference_order", args.data_root/"lib/order.py")
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    start = time.perf_counter()
    sources = prepare(args.data_root, args.output, args.frames, args.shooting_frames, not args.no_quenched)
    if args.only:
        allowed = set(args.only.split(","))
        sources = [s for s in sources if s["name"] in allowed]
    if args.limit_sources:
        sources = sources[:args.limit_sources]
    if not sources:
        raise ValueError("No selected trajectory sources")
    manifest = {"source_project": str(args.data_root), "cutoff_A": 5., "first_shell_cutoff_A": 3.5,
                "n_jobs": args.n_jobs, "sampling": "Evenly spaced stored-frame indices including endpoints",
                "sources": sources, "snapshot_checksums": [],
                "reference_order_sha256": hashlib.sha256((args.data_root/"lib/order.py").read_bytes()).hexdigest(),
                "versions": {p: version(p) for p in ("SolOrder", "numpy", "scipy", "MDAnalysis", "numba", "joblib", "freud-analysis", "h5py")}}
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    rows, groups, checks = [], [], []
    for source in sources:
        print(f"SOURCE {source['name']}: {len(source['selected_frames'])}/{source['frames']} frames", flush=True)
        u = snapshot_template(source["topology"]) if source["full_atom"] else oxygen_universe(source["atoms"])
        for number, (index, positions, dimensions, time_ns) in enumerate(iter_frames(source)):
            tick = time.perf_counter()
            u.atoms.positions = positions
            u.dimensions = dimensions
            op = OrderParameter(u, cutoff=5, selection="name O OW", hydrogen_selection="name H H1 H2 HW1 HW2", n_jobs=args.n_jobs)
            np.testing.assert_array_equal(op.pos, np.asarray(positions, dtype=np.float32)[op.atom_indices])
            np.testing.assert_array_equal(op.dimensions, np.asarray(dimensions, dtype=np.float32))
            values, pair_angles = fields(op, source["full_atom"])
            np.testing.assert_array_equal(u.atoms.positions, np.asarray(positions, dtype=np.float32))
            label_result = op.chill_plus()
            original_labels, _ = reference.chill_plus(op.pos/10, op.dimensions[:3]/10)
            mismatches = int(np.count_nonzero(label_result.labels != original_labels))
            cn_reference = np.bincount(reference.neighbor_pairs(op.pos/10, op.dimensions[:3]/10)[0].ravel(), minlength=op.n)
            core = largest_core(op, label_result.labels)
            all_core = largest_core(op, label_result.labels, (1, 2, 3))
            check = dict(source=source["name"], frame=index, centers=op.n,
                         snapshot_coordinate_match=True, label_mismatches=mismatches,
                         cn_mismatches=int(np.count_nonzero(values["cn_3p5"] != cn_reference)),
                         voronoi_relative_volume_error=abs(values["voronoi_volume"].sum()/op.volume-1),
                         hb_mismatches=None, parallel_max_abs_error=0.,
                         bound_violations=bound_violations(values),
                         inf_count=sum(int(np.isinf(a).sum()) for a in values.values()),
                         nan_count=sum(int(np.isnan(a).sum()) for a in values.values()))
            if source["full_atom"]:
                check["hb_mismatches"] = int(np.count_nonzero(values["HB"] != reference_hb(op, reference)))
            if number == 0:
                serial = OrderParameter(u, cutoff=5, selection=op.selection, hydrogen_selection=op.hydrogen_selection)
                serial_values, serial_angles = fields(serial, source["full_atom"])
                for name, value in values.items():
                    other = serial_values[name]
                    if not np.array_equal(np.isnan(value), np.isnan(other)):
                        raise AssertionError(f"Serial/parallel NaN mismatch: {source['name']}/{name}")
                    finite = np.isfinite(value) & np.isfinite(other)
                    error = float(np.max(np.abs(value[finite]-other[finite]))) if finite.any() else 0.
                    check["parallel_max_abs_error"] = max(check["parallel_max_abs_error"], error)
                    np.testing.assert_allclose(value, other, rtol=1e-10, atol=1e-10, equal_nan=True)
                np.testing.assert_array_equal(label_result.labels, serial.chill_plus().labels)
            else:
                check["parallel_max_abs_error"] = None
            checks.append(check)
            row = dict(source=source["name"], family=source["family"], frame=index,
                       stored_frames=source["frames"], progress=index/max(1, source["frames"]-1), time_ns=time_ns,
                       waters=op.n, full_atom=source["full_atom"], box_approximate=source["box_approximate"],
                       density_Aminus3=op.density(), box_x_A=float(dimensions[0]),
                       largest_Ih_Ic_core=int(core.sum()), largest_ice_including_interface=int(all_core.sum()))
            for code, label in enumerate(LABELS):
                row["chill_plus_"+label] = int(np.count_nonzero(label_result.labels == code))
            classic = op.chill()
            for code, label in enumerate(LABELS):
                row["chill_"+label] = int(np.count_nonzero(classic.labels == code))
            for name, array in values.items():
                row[f"{name}_core_auc"] = core_auc(array, core)
                for statistic, value in summarize(array).items():
                    row[f"{name}_{statistic}"] = value
            if pair_angles is not None:
                for statistic, value in summarize(pair_angles.angles).items():
                    row[f"Ptheta_{statistic}"] = value
                check["Ptheta_integral"] = float(np.nansum(pair_angles.probability_density*np.diff(pair_angles.bin_edges)))
            rows.append(row)
            for group_name, mask in (("largest_core", core), ("outside_core", ~core), ("all_Ih_Ic", np.isin(label_result.labels, (1, 2)))):
                for name, array in values.items():
                    record = dict(source=source["name"], frame=index, progress=row["progress"], group=group_name,
                                  group_size=int(mask.sum()), parameter=name, **summarize(array[mask]))
                    groups.append(record)
            manifest["snapshot_checksums"].append(dict(source=source["name"], frame=index,
                analysis_oxygen_sha256=hashlib.sha256(op.pos.tobytes()+op.dimensions.tobytes()).hexdigest(),
                sha256=hashlib.sha256(np.asarray(positions, dtype="<f8").tobytes()+np.asarray(dimensions, dtype="<f8").tobytes()).hexdigest()))
            print(f"  {index:5d} core={core.sum():4d} Q={row['Q_mean']:.4f} label_delta={mismatches} "
                  f"elapsed={time.perf_counter()-tick:.2f}s", flush=True)
        write_csv(args.output/"frames.csv", rows)
        write_csv(args.output/"phase_conditioned.csv", groups)
        write_csv(args.output/"checks.csv", checks)
        (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    manifest["elapsed_seconds"] = time.perf_counter()-start
    manifest["analyzed_frames"] = len(rows)
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    failures = [c for c in checks if c["label_mismatches"] or c["cn_mismatches"] or c["hb_mismatches"] or c["bound_violations"] or c["inf_count"] or c["voronoi_relative_volume_error"] > 1e-5]
    manifest["failed_frames"] = len(failures)
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"Finished {len(rows)} sampled frames in {manifest['elapsed_seconds']:.1f}s; {len(failures)} frames with failed checks", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
