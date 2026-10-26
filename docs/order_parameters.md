# Parallel structural order parameters

SolOrder analyzes a copied MDAnalysis trajectory frame. Per-oxygen result arrays
follow `analysis.atom_indices`, which contains the original topology indices.
Atom IDs and residue IDs need not be consecutive or unique. Water calculations
require one selected oxygen and two selected hydrogens per residue.

```python
import MDAnalysis as mda
from SolOrder import OrderParameter, eFingerprint

u = mda.Universe("water.gro", "water.xtc")
analysis = OrderParameter(
    u, frame=0, cutoff=5.0, n_jobs=4, backend="threading",
    selection="name OW", hydrogen_selection="name HW1 HW2",
)

number_density = analysis.density()             # selected oxygens / Angstrom^3
cn = analysis.coordination_number(cutoff=3.5)
smooth_cn = analysis.smooth_coordination(r0=3.5, cutoff=5.0)
d5 = analysis.d5()                               # Angstrom
q_tetrahedral = analysis.OTO()
q_tetrahedral_averaged = analysis.tetrahedral_order(averaged=True)
s_tetrahedral = analysis.TTO()
lsi = analysis.LSI()                             # Angstrom^2
f4 = analysis.F4(cutoff=3.5)

q4 = analysis.steinhardt(l=4)
q6 = analysis.steinhardt(l=6)
q3_averaged = analysis.steinhardt(l=3, averaged=True)
w6 = analysis.steinhardt(l=6, invariant="w")
w6_normalized = analysis.steinhardt(l=6, invariant="w", normalized=True)
connections = analysis.crystalline_connections(l=6, threshold=0.5)

hb = analysis.hydrogen_bonds(cutoff=3.5, angle=30.0)
# hb.donated, hb.accepted, hb.total: one count per selected oxygen.
# hb.bonds columns: donor oxygen row, hydrogen topology index, acceptor oxygen row.

entropy = eFingerprint(u, frame=0, cutoff=5.0, local=False, n_jobs=4)
e = entropy.Entropy()                            # units of k_B
le = entropy.localEntropy()

# Advancing the live Universe does not alter the analysis snapshot.
analysis.refresh(frame=10)
q6_frame10 = analysis.steinhardt(l=6)
```

`StructuralOrder` exposes the same new lowercase methods without the historical
`OTO`, `TTO`, `F4`, `LSI` and `MinimumAngle` wrappers. The wrappers continue to
set `tetra_orient`, `tetra_trans`, `f4`, `lsi` and `minAngles`, respectively, and
now also return their arrays. `get_dihedrals()` stores unique-pair dihedrals in
`dihedrals`. `refresh()` invalidates caches and previously stored results.

## Parallel execution

`n_jobs=1` executes serially, positive values request that many workers, and
`n_jobs=-1` requests all available workers. `backend="threading"` shares the
snapshot arrays. `backend="loky"` uses joblib processes. Worker counts are capped
by the number of input items. Kernel outputs retain selection order.

The helper in `SolOrder.parallel` splits atoms or neighbor edges into chunks,
runs the first chunk on the caller to compile Numba kernels, then dispatches the
remaining chunks through joblib. Kernels use `njit(cache=True, nogil=True)` and
no internal parallel thread pool. Neighbor geometry uses MDAnalysis's compiled
spatial search and minimum-image functions. Spherical harmonics use SciPy;
Numba performs the reductions and Wigner contractions. Scalar density is a
simple volume division and does not need parallel execution.

First calls include JIT compilation and neighbor-list construction. Small
systems may run faster with `n_jobs=1`; parallelism is not a guaranteed speedup.
Do not access or advance a Universe inside user-supplied parallel work.

## Definitions and boundaries

- A valid six-component periodic box is required. Orthorhombic and triclinic
  geometry, including coordinates shifted by multiple cells, are supported.
- Density is global **number density** of the selected particles, not mass
  density or Voronoi density.
- Cutoff coordination counts other selected atoms with `distance <= cutoff`.
  Smooth coordination uses exponents 12 and 28 by default, with the analytic
  limit `12/28` at `distance == r0`. Its infinite tail is truncated at the search
  cutoff; increase that cutoff to check convergence.
- Nearest-neighbor quantities search within the configured cutoff. `d5()` gives
  NaN when fewer than five neighbors are found; increase the cutoff to obtain
  the actual fifth neighbor. Tetrahedral order needs four neighbors. LSI needs
  at least one neighbor inside 3.7 Angstrom and the first neighbor outside it.
  Missing information produces NaN, rather than a plausible finite value.
- LSI includes the gap from the last inner-shell neighbor to the first outer
  neighbor. Translational tetrahedral order uses **squared** radial deviations.
- Local averages include the central atom and use uniform neighbor weights.
  NaNs propagate. Averaged Steinhardt order averages the complex `q_lm` vectors
  before forming the invariant; it does not average scalar `q_l` values.
- `q_l` supports integer `0 <= l <= 32`; the factorial-sum `w_l` implementation
  is limited to `l <= 12`. Normalization applies only to `w_l`. Coincident atoms
  have undefined bond orientations and cause a clear error. Crystalline
  connection thresholds must be chosen for the system; the count alone does
  not assign an ice phase.
- Hydrogen bonds use O--O distance at most 3.5 Angstrom and donor H--O--O angle
  strictly below 30 degrees by default. This is the angle at the donor oxygen,
  not the complementary O--H--O angle. Counts represent donor-hydrogen-acceptor
  bonds; each contributes to one donated and one accepted count. Hydrogen
  association uses residue membership and the explicit hydrogen selection.
- F4 independently selects each molecule's hydrogen farthest from the opposite
  oxygen, using minimum-image O--O and O--H
  vectors, and a per-call cutoff override. That override does not alter other
  calculations' cutoffs. `MinimumAngle()` reports each oxygen's minimum angle
  with its own nearest oxygen; it no longer copies one molecule's result onto
  a neighbor with a different nearest-neighbor environment.
- Entropy honors `frame` and `local`. `include_hydrogen=True` uses all atoms for
  radial density while returning and averaging only oxygen-center values.
  The integration sphere must fit within half the smallest periodic cell
  height. An isolated center has undefined local normalization (NaN); global
  normalization uses the finite zero-neighbor limit. Default frame is now 0.

## Validation and references

Tests cover known FCC q4/q6 and normalized w4/w6 values, the Legendre addition
formula as an independent check of q_l, rotations, periodic translations,
triclinic boxes, direct hydrogen-bond geometry, independent pair-entropy
integration, and agreement between serial, thread and process execution.
The [GenIce regression report](validation/genice.md) adds checks on eleven
ice/clathrate structures. These tests do not establish phase-classification
accuracy on thermally sampled trajectories.

Implementation conventions follow the manuscript and these primary references:

- [Steinhardt and averaged order definitions](https://freud.readthedocs.io/en/latest/modules/order.html)
- [SciPy spherical-harmonic angle conventions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html)
- [MDAnalysis periodic geometry](https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html)

Periodic Voronoi geometry, CHILL/CHILL+, the remaining geometric manuscript
parameters and non-neural descriptor APIs are documented in
[extended parameters](extended_parameters.md). Complex learned models are
deferred at the user's request.
