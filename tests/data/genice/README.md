# GenIce regression structures

These are raw Gromacs outputs generated with **GenIce2 2.2.13.3**, TIP3P water,
seed 12345, and `PYTHONHASHSEED=0`. Atom names are the original `O` and `H`.
Gromacs stores coordinates in nm; MDAnalysis reads them in Angstrom.

Included lattices: Ih, Ic, II, III, V, VI, VII, VIII, XII, XVII, and the empty sI
clathrate water framework (no methane or other guest). They are as-generated
structures, without energy minimization, thermalization or MD sampling.

`manifest.json` records each exact command, repeat count, dimensions, water
count, file SHA-256, and the generator environment's dependency versions.
Repetition is increased if necessary so the shortest cell height is >=12 A.

Recreate them from the repository root using two separate environments:

```sh
.conda-env/bin/python -m venv .genice-env
.genice-env/bin/python -m pip install genice2==2.2.13.3
.conda-env/bin/python scripts/generate_genice_fixtures.py
.conda-env/bin/python -m pytest -q tests/test_genice.py
.conda-env/bin/python scripts/validate_genice.py
```

GenIce2 requires NumPy 2; SolOrder currently uses NumPy 1. Keeping the generator
separate avoids changing the analysis environment. GenIce is not a runtime or
test dependency: the checked-in structures support offline tests.

Source: https://github.com/genice-dev/GenIce2 (MIT-licensed generator).
