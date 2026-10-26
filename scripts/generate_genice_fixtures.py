"""Generate reproducible GenIce2 fixtures without importing it into SolOrder.

Run with SolOrder's Python and a separate GenIce2 executable. Coordinates remain
in the generator's original Gromacs format (O/H names, nm in file).
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.mdamath import triclinic_vectors


PHASES = {
    "Ih": "1h", "Ic": "1c", "II": "2", "III": "3", "V": "5", "VI": "6",
    "VII": "7", "VIII": "8", "XII": "12", "XVII": "17", "sI": "sI",
}
# Benchmark-only lattices (docs/benchmark_design.md, dataset D1); not test fixtures.
EXTRA_PHASES = {
    "IV": "4", "IX": "9", "XI": "11", "XIII": "13", "XIV": "14", "XVI": "16",
    "sII": "sII", "sH": "sH",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genice", type=Path, default=root/".genice-env/bin/genice2")
    parser.add_argument("--output", type=Path, default=root/"tests/data/genice")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--phases", help="Comma-separated phase names (default: the test fixtures)")
    args = parser.parse_args()
    known = {**PHASES, **EXTRA_PHASES}
    phases = args.phases.split(",") if args.phases else list(PHASES)
    unknown = [p for p in phases if p not in known]
    if unknown:
        raise SystemExit(f"Unknown phases: {unknown}; choose from {list(known)}")
    args.output.mkdir(parents=True, exist_ok=True)
    generator = args.genice.resolve()
    env = dict(os.environ, PYTHONHASHSEED="0")
    version = subprocess.check_output([str(generator), "--version"], text=True, env=env).strip()
    manifest = {"generator": version, "seed": args.seed, "pythonhashseed": 0,
                "format": "gromacs", "water": "tip3p", "selection": "name O",
                "hydrogen_selection": "name H", "structures": {}}
    for phase in phases:
        lattice = known[phase]
        target = args.output/f"{phase}.gro"
        # Keep the entropy sphere (5 Angstrom) well within every periodic height.
        for repeat in range(3, 8):
            command = [str(generator), lattice, "--rep", *([str(repeat)]*3),
                       "--seed", str(args.seed), "--water", "tip3p", "--format", "gromacs", "--quiet"]
            result = subprocess.run(command, text=True, capture_output=True, env=env, check=True)
            target.write_text(result.stdout)
            universe = mda.Universe(str(target))
            cell = triclinic_vectors(universe.dimensions, dtype=np.float64)
            volume = abs(np.linalg.det(cell))
            heights = volume / np.linalg.norm(np.cross(cell[[1, 2, 0]], cell[[2, 0, 1]]), axis=1)
            if min(heights) >= 12:
                break
        else:
            raise RuntimeError(f"Could not create a sufficiently large box for {phase}")
        manifest["structures"][phase] = {
            "lattice": lattice, "file": target.name, "rep": [repeat]*3,
            "command": ["genice2", *command[1:]],
            "molecules": len(universe.select_atoms("name O")),
            "dimensions_angstrom_degrees": universe.dimensions.tolist(),
            "minimum_cell_height_angstrom": float(min(heights)),
            "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        }
        print(f"{phase}: {manifest['structures'][phase]['molecules']} waters; min height {min(heights):.2f} A", flush=True)
    manifest["generator_dependencies"] = subprocess.check_output(
        [str(generator.parent/"python"), "-m", "pip", "freeze"], text=True, env=env
    ).splitlines()
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    main()
