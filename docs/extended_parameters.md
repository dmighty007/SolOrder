# Extended parameters and descriptor APIs

All methods below operate on the same snapshot as `OrderParameter` or
`StructuralOrder`. Numeric distances use Angstrom, and per-center outputs follow
`op.atom_indices`. Advance a trajectory with `op.refresh(frame)`.

```python
from SolOrder import OrderParameter

op = OrderParameter(u, cutoff=5., n_jobs=4)
volume = op.voronoi_geometry().volume
local_density = op.voronoi_density()
asphericity = op.voronoi_asphericity()

angles = op.water_angular_distribution(cutoff=3.5, bins=90)
# angles.angles and angles.pairs: one entry per unique oxygen pair.
# angles.bin_edges, angles.counts, angles.probability_density: radians/radian^-1.

csp = op.centrosymmetry(n_neighbors=4)
csp_averaged = op.centrosymmetry(4, averaged=True)
distance_op = op.neighbor_distance_parameter(f_alpha=1, f_beta=2, f_gamma=1)
fourier_op = op.angular_fourier(a=3., radial=0)
bond_angle = op.bond_angle_parameter(n1=2, n2=1, phi=0.)
angle_fft = op.angle_histogram_parameter(frequencies=[1, 2, 3], bins=180)
a_parameters = op.neighborhood_parameters(m=1)  # N x 3: A1, A2, A3.

chill = op.chill()
chill_plus = op.chill_plus()
print(chill_plus.names)
# Result also includes labels, staggered/eclipsed counts, correlations,
# directed neighbor offsets and indices.
```

## Periodic Voronoi measures

freud/Voro++ computes periodic cells from all selected atoms, including a full
triclinic cell. Volumes sum to the simulation volume. Results are cell volume,
face area, inverse volume and area^3/(36 pi volume^2). They are cached until
refresh; coincident positions are rejected. Voronoi uses its compiled backend,
which controls its own workers. It does not launch an outer joblib pool.

## Angle distribution

The angle for each unique O--O pair is the smallest of the four donor-oxygen
H--O--O angles, as in the manuscript. All pairs within the explicit cutoff are
included. The probability density integrates to one over [0, pi]. Empty samples
return zero counts and NaN density. This differs from `MinimumAngle()`, which
uses only each oxygen's own nearest partner.

## Geometric conventions

- Conventional centrosymmetry sums the N/2 smallest squared pair-vector sums,
  as in LAMMPS. `pairing="rank"` instead pairs distance-ranked first-half and
  second-half neighbors, matching the manuscript's expression. N must be even;
  insufficient neighbors give NaN. These two definitions can differ in distorted
  environments.
- Neighbor-distance parameters average a product of three radial functions over
  distinct pairs of neighbors. Integer arguments specify powers of distance;
  callables receive distance arrays. Alpha/beta factors are symmetrized by
  default; `symmetrize=False` uses distance-ranked pair ordering. The opposite
  neighbor-neighbor distance uses minimum-image geometry too.
- Angular Fourier `weighting="min"` follows the manuscript's minimum-radius
  radial factor. `weighting="product", radial_beta=...` selects the product of
  separate radial factors used in MALIO; `symmetrize=False` preserves its rank
  ordering. The default zero radial power gives unit weights.
- Bond-angle B averages cos(n2 theta + phi)^n1. This fills the manuscript's
  omitted function definition using the published source.
- Histogram H uses equal bins on [0, pi], pair-normalized counts and DFT
  amplitude channels. `power=True` returns squared amplitudes as in MALIO.
  Frequencies are integer rFFT indices up to bins/2. Averaging occurs after the
  amplitude/power calculation. The exact pi endpoint belongs to the last bin.
- A1/A2/A3 use the author's MALIO convention. For each i--j bond, other
  neighbors k of i are ranked by r_ik+r_jk; up to m are used. A1 and A2 average
  squared vector sums over j, whereas A3 squares the combined sum over j.
  This resolves the duplicated A2/A3 manuscript equations. This ranking does
  not impose a second cutoff-neighbor intersection. Outputs have units A^2.
- Neighbor-pair statistics require at least two neighbors; missing environments
  give NaN. `averaged=True` includes the central atom and uses uniform neighbor
  weights. NaNs propagate.

## CHILL and CHILL+

CHILL uses four closest oxygens inside the explicit search cutoff, normalized
q3 correlations, a staggered bound below -0.8 and eclipsed interval
(-0.2,-0.05). Interfacial classification uses the original neighbor-dependent
rules. Increase the search cutoff if four neighbors are not found.

CHILL+ follows OVITO's cutoff convention. All neighbors inside the cutoff enter
q3m. Every neighbor with nonzero harmonic norm can contribute a correlation,
including undercoordinated neighbors. A center receives a structure label only
when it has exactly four neighbors. Its eclipsed interval is (-0.35,0.25).
This matters at surfaces and in disordered or overcoordinated environments.

Codes are 0 other, 1 hexagonal ice, 2 cubic ice, 3 interfacial ice, 4 hydrate,
and 5 interfacial hydrate. Code 0 includes structures outside these rules; it
is not a general liquid label. CHILL does not assign hydrate codes. Classifiers
need only oxygen positions. Thresholds are fixed algorithm definitions, with
no learned weights or fitted phase probabilities.

## Model-dependent molecular quadrupole

```python
qt = op.quadrupole_parameter(dipole=2.35, negative_site_distance=0.15)
```

Supply dipole magnitudes in the units used by your water model and d_OM in
Angstrom. The three-charge-model expression is
Q_T = 3 (r_OH sin(theta))^2 mu / [4 (r_OH cos(theta) - d_OM)], with theta half
the H-O-H angle. Output units are dipole units times Angstrom. The mean of the
two O-H lengths approximates flexible geometries. Both inputs accept one value
per water. The default d_OM=0 reproduces the manuscript expression; model-specific
charge placement requires its actual offset. This molecular moment does not
classify the oxygen lattice.

## Optional descriptors

```sh
python -m pip install -e ".[descriptors]"
```

The optional dependency range uses ASE <3.27 because DScribe 2.1's legacy
constructor call is incompatible with newer ASE constructors. Tests use ASE
3.26.0 and DScribe 2.1.2.

```python
soap = op.soap(symbols="O", r_cut=5., n_max=8, l_max=6, sigma=.5)
acsf = op.atomic_descriptor("acsf", symbols="O", r_cut=5.,
                            g2_params=[[1., 0.], [.5, 2.]])
lmbtr = op.atomic_descriptor(
    "lmbtr", symbols="O", geometry={"function": "distance"},
    grid={"min": 0, "max": 6, "sigma": .1, "n": 60},
    weighting={"function": "exp", "scale": 1., "threshold": 1e-3},
)
mbtr = op.atomic_descriptor(
    "mbtr", symbols="O", geometry={"function": "distance"},
    grid={"min": 0, "max": 6, "sigma": .1, "n": 60},
    weighting={"function": "exp", "scale": 1., "threshold": 1e-3},
)
cm = op.atomic_descriptor("coulomb", symbols="O", n_atoms_max=op.n)
sine = op.atomic_descriptor("sine", symbols="O", n_atoms_max=op.n)
ewald = op.atomic_descriptor("ewald", symbols="O", n_atoms_max=op.n)
graph = op.graph_descriptor(cutoff=3.5)
```

Descriptors include the selected atoms only. Supply chemical symbols explicitly
for mixed selections; atom names are not guessed. SOAP/LMBTR/ACSF return per-center
rows and support `centers=[...]` and joblib chunks; each chunk has its own
backend object. Global MBTR and matrix descriptors return one flattened vector
from the compiled backend. Matrices use atomic numbers as defined by DScribe,
not water partial charges. Coulomb matrices describe a nonperiodic configuration;
sine and Ewald matrices include the supplied periodic cell. Matrix permutation,
padding and descriptor settings are passed explicitly to DScribe.

The graph result exposes directed minimum-image edges, topology-index mapping,
node degree/mean-distance/distance-standard-deviation features and edge distances.
These scalar features are invariant to rigid rotation. No graph neural network
is invoked.

## Optional template matching

```sh
python -m pip install -e ".[ptm]"
```

```python
result = op.template_matching(
    templates=("CUBIC_DIAMOND", "HEX_DIAMOND"), rmsd_cutoff=.1,
)
```

OVITO's compiled PTM backend supplies structure codes, type names and RMSD.
Default diamond templates target ice Ic/Ih oxygen networks. Other OVITO templates
can be selected explicitly. OVITO manages its own worker pool; `op.n_jobs` does
not configure that backend. This does not identify proton ordering or supply
all ice-polymorph templates.

DeepIce, GCIceNet and the SOAP autoencoder remain deferred at the user's request.

## References

- [Takahashi (2023): geometric parameter definitions](https://doi.org/10.1039/D2CP03696G)
- [Author's MALIO source: A/P/N, radial and histogram conventions](https://github.com/Kazuaki-Z-Takahashi/MALIO)
- [OVITO CHILL+ definitions](https://www.ovito.org/manual/reference/pipelines/modifiers/chill_plus.html)
- [OVITO CHILL+ reference implementation](https://gitlab.com/stuko/ovito/-/blob/master/src/ovito/particles/modifier/analysis/chill_plus/ChillPlusModifier.cpp)
- [Centrosymmetry conventions](https://www.ovito.org/manual/reference/pipelines/modifiers/centrosymmetry.html)
- [freud periodic Voronoi API](https://freud.readthedocs.io/en/latest/modules/locality.html#freud.locality.Voronoi)
- [DScribe descriptor definitions](https://singroup.github.io/dscribe/latest/)
