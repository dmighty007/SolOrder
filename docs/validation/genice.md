# Validation against GenIce structures

Implementation coverage: **25 implemented** in the original
25-calculation scope; complex learned models are deferred. See the [order-parameter status table](../order_parameter_status.md)
for APIs, validation coverage and the remaining manuscript methods.

The regression set contains **11 structures and 6,696 water molecules**,
generated with GenIce2 2.2.13.3, TIP3P, seed 12345 and PYTHONHASHSEED=0.
It spans orthorhombic and triclinic cells: Ih, Ic, II, III, V, VI, VII, VIII,
XII, XVII, and an empty sI clathrate water framework. Original O/H names and
Gromacs coordinates are preserved. Each periodic cell is large enough for the
5-Angstrom entropy integration sphere.

These are static, as-generated structures. Values are fixture-specific
regression results, not thermally equilibrated phase averages or universal
classification thresholds. All reported means are arithmetic means over oxygen
centers; F4 is first averaged over each center's neighbors.

| Structure | Waters | Mean Q | Mean CN (3.5 A) | Mean F4 (3.5 A) | Mean total H-bonds |
|---|---:|---:|---:|---:|---:|
| Ih | 432 | 0.998073 | 4.0000 | -0.493395 | 4 |
| Ic | 216 | 0.998173 | 4.0000 | -0.987278 | 4 |
| II | 324 | 0.827606 | 5.5000 | -0.167193 | 4 |
| III | 324 | 0.857022 | 4.5432 | -0.157922 | 4 |
| V | 756 | 0.798020 | 5.8228 | -0.110373 | 4 |
| VI | 270 | 0.687173 | 11.4593 | -0.004466 | 4 |
| VII | 432 | 0.403271 | 14.0000 | -0.358224 | 4 |
| VIII | 1728 | 0.332222 | 13.0000 | -0.476305 | 4 |
| XII | 648 | 0.787493 | 9.5586 | -0.072985 | 4 |
| XVII | 324 | 0.971669 | 4.0000 | 0.497846 | 4 |
| sI | 1242 | 0.982343 | 4.0000 | 0.908459 | 4 |

Full means and standard deviations, including smooth coordination, d5, LSI,
q3/q4/q6, averaged q3/q4, normalized w6, crystalline connections, and local/global
entropy, are in [genice_summary.csv](genice_summary.csv).

## Checks

The GenIce tests check:

- Fixture SHA-256 and molecule counts against the generation manifest.
- Exactly two donated and two accepted hydrogen bonds at every water in all
  eleven structures, with no duplicated donor-hydrogen-acceptor entries.
- Coordination, d5, density and tetrahedral formulas against direct reference
  calculations. Periodic distances are independently checked by enumerating 125
  lattice images for twelve distributed center atoms in each structure.
- q3/q4/q6 against the Legendre addition formula and LSI against sorted shell
  gaps. F4 is checked against independent outer-hydrogen selection and
  MDAnalysis dihedrals for all oxygen centers.
- Finite outputs, serial/thread agreement and whole-water periodic translation
  invariance. Process-backend agreement is additionally checked on Ih and
  monoclinic V. Local and global entropy are checked in serial and parallel.
- The synthetic/unit suite separately checks independent entropy integration,
  FCC q4/q6/w4/w6 values, rotations, no-neighbor behavior, and other edge cases.

### F4 discrepancy found and corrected

The inherited algorithm selected hydrogens using the global minimum and maximum
of four cross-water distances. If both extrema belonged to one molecule, the
other molecule's selected hydrogen could be wrong. GenIce II and III exposed
per-center errors up to about 0.10 and 0.33, respectively; VIII showed errors up
to about 0.44. The kernel and compatibility helper now choose the hydrogen
farthest from the opposite oxygen **independently for each molecule**. The
independent dihedral reference passes across all eleven structures.

This follows the outermost-hydrogen definition of F4 in
[Jing et al. (2024)](https://doi.org/10.1186/s13065-024-01191-6).
It changes previously calculated F4 values in affected environments.

### Neighbor-shell ambiguity

Ice VII and VIII can have degenerate fourth/fifth-neighbor distances because
multiple networks interpenetrate. A four-neighbor tetrahedral score may depend
on which equally distant neighbors are selected. SolOrder sorts by distance,
then atom index. Tests compare exact four-neighbor identities only where the
fourth neighbor is unambiguous; cutoff-based bond-order tests include every
neighbor within the cutoff. A low Q score here does not by itself mean that
hydrogen bonds or the ice structure are invalid. Every structure in this set
still has exactly four geometric hydrogen bonds per water.

## Reproduction

The normal test suite uses bundled structures and requires no GenIce install:

```sh
.conda-env/bin/python -m pytest -q
.conda-env/bin/python scripts/validate_genice.py
```

For generation commands, environment versions and input hashes, see
[the fixture manifest](../../tests/data/genice/manifest.json) and
[fixture instructions](../../tests/data/genice/README.md). The generation script
has also been rerun into a separate directory to check byte-for-byte input
reproducibility in the recorded environment.

Generator: [GenIce2](https://github.com/genice-dev/GenIce2).

## Extended-parameter validation

The suite now passes **136 tests**, covering the original methods, remaining
geometric parameters, and optional descriptor/template integrations. New tests
include exact cubic Voronoi volumes/areas/asphericity, triclinic translation,
centrosymmetry pair-weight sums, independent angular/triangle formulas, normalized
histograms, A1/A2/A3 vector sums and model-dependent quadrupole geometry.

CHILL correctly identifies the Ih and Ic fixtures. CHILL+ matches OVITO's
reference labels on every GenIce input below, and on three additional perturbed
Ih slabs with undercoordinated surfaces. Its local q3m uses all cutoff neighbors,
including vectors at undercoordinated neighboring sites; only centers with
exactly four neighbors receive a phase label. Labels describe local motifs: some
III and XVII centers receive interface labels, as shown below. These outputs do
not constitute a complete ice-polymorph classifier.

| Structure | Mean Voronoi volume (A^3) | Mean asphericity | CHILL+ labels and counts |
|---|---:|---:|---|
| Ih | 32.4896 | 2.2039 | Ih: 432 |
| Ic | 32.4896 | 2.1605 | Ic: 216 |
| II | 25.3542 | 1.6488 | other: 324 |
| III | 25.6570 | 1.6217 | other: 216; ice_interface: 108 |
| V | 24.1084 | 1.5032 | other: 756 |
| VI | 21.7701 | 1.3950 | other: 270 |
| VII | 18.6815 | 1.3284 | other: 432 |
| VIII | 18.3588 | 1.3319 | other: 1728 |
| XII | 20.7615 | 1.4586 | other: 648 |
| XVII | 33.9664 | 2.3767 | other: 266; hydrate_interface: 58 |
| sI | 37.5980 | 2.6988 | hydrate: 1242 |

The CSV now also contains Voronoi geometry, centrosymmetry, neighbor-distance,
angular Fourier, B, A1/A2/A3, angle-histogram channels, water-pair angles and both
CHILL classification counts. Optional tests exercise SOAP, MBTR/LMBTR, ACSF,
Coulomb/sine/Ewald descriptors, graph features and diamond-template PTM.

Reference-backend versions used here: freud 3.4.0, OVITO 3.15.5, DScribe 2.1.2,
ASE 3.26.0. DScribe/PTM tests skip when their extras are absent; they ran in this
validation environment. Complicated learned models remain deferred.

See [extended API conventions](../extended_parameters.md) for definitions,
radial-weight and pairing choices, species requirements and optional installs.
