# Order-parameter implementation status

All **25 calculations in the original scope are implemented**. Additional
non-neural manuscript methods are available below, including optional descriptor
and template-matching backends. DeepIce, GCIceNet and the SOAP autoencoder are
deferred at the user's request.

`op` denotes `OrderParameter(...)`; its lowercase structural methods also exist
on `StructuralOrder(...)`. `entropy` denotes `eFingerprint(...)`. Implementation
means an executable package API exists. Validation describes the checks actually
performed, rather than establishing phase-classification accuracy.

## Original 25-calculation scope

| # | Order parameter / calculation | Status | API or remaining work | Validation |
|---:|---|---|---|---|
| 1 | Global number density | Implemented | `op.density()`; selected atoms per A^3, not mass density | Direct box-volume reference; GenIce |
| 2 | Cutoff coordination number | Implemented | `op.coordination_number(cutoff=3.5)` | Direct distance counts; GenIce |
| 3 | Smooth coordination number | Implemented | `op.smooth_coordination(r0=3.5, n=12, m=28)` | Direct switching-function reference and analytic limit; GenIce output summary |
| 4 | Orientational tetrahedral order Q | Implemented | `op.tetrahedral_order()` or `op.OTO()` | Ideal tetrahedron and direct formulas; GenIce, excluding ambiguous fourth-neighbor identities |
| 5 | Locally averaged tetrahedral order | Implemented | `op.tetrahedral_order(averaged=True)` | Unit test of neighbor averaging |
| 6 | Translational tetrahedral order | Implemented | `op.translational_order()` or `op.TTO()` | Ideal tetrahedron and direct radial-variance reference; GenIce |
| 7 | Local pair-entropy fingerprint | Implemented | `entropy.Entropy()`; local/global density normalization | Independent integration reference; GenIce finite-output and serial/parallel checks |
| 8 | Locally averaged entropy fingerprint | Implemented | `entropy.localEntropy()`; uniform oxygen-neighbor weights | Neighbor-average unit tests; GenIce serial/parallel checks |
| 9 | Local Structure Index (LSI) | Implemented | `op.local_structure_index()` or `op.LSI()` | Direct shell-gap variance; GenIce |
| 10 | Four-body order F4 | Implemented | `op.F4(cutoff=3.5)` | Independent outer-hydrogen and dihedral reference; GenIce; hydrogen-selection bug corrected |
| 11 | Fifth-nearest oxygen distance d5 | Implemented | `op.d5()`; NaN when search cutoff contains fewer than five neighbors | Sorted-distance reference; GenIce |
| 12 | Voronoi number density | Implemented | `op.voronoi_density()`; periodic freud/Voro++ cells | Exact cubic cells; volume conservation on all GenIce phases; triclinic translation |
| 13 | Voronoi asphericity | Implemented | `op.voronoi_asphericity()`; area^3/(36 pi volume^2) | Exact cube area and asphericity; finite physical results on all GenIce phases |
| 14 | Hydrogen-bond count | Implemented | `op.hydrogen_bonds()` returns donated, accepted and total counts | Direct geometry; every GenIce water has two donated and two accepted bonds |
| 15 | Water-water angular distribution P(theta) | Implemented | `op.water_angular_distribution(cutoff=3.5, bins=90)`; angles, unique pairs, counts and normalized density | Independent H-O-O geometry; histogram normalization; GenIce serial/parallel checks |
| 16 | Steinhardt q_l | Implemented | `op.steinhardt(l=6)` | FCC values, Legendre reference and rotations; GenIce q3/q4/q6 |
| 17 | Averaged Steinhardt q_l | Implemented | `op.steinhardt(l=6, averaged=True)` | Complex-vector averaging and rotation tests; GenIce serial/parallel checks |
| 18 | Third-order Steinhardt w_l | Implemented | `op.steinhardt(l=6, invariant="w")`; l <= 12 | Shared Wigner contraction exercised by normalized-w tests; no separate raw-w reference test |
| 19 | Normalized Steinhardt w_l | Implemented | `op.steinhardt(l=6, invariant="w", normalized=True)` | FCC w4/w6, rotations; GenIce normalized-w6 serial/parallel checks |
| 20 | Crystalline connection count | Implemented | `op.crystalline_connections(l=6, threshold=0.5)`; system-specific threshold | Unit count check; GenIce output summary; phase thresholds not calibrated |
| 21 | Centrosymmetry parameter | Implemented | `op.centrosymmetry(n_neighbors=12, pairing="conventional")`; `pairing="rank"` for manuscript rank pairing; optional averaging | Exact centrosymmetric and distorted environments; direct pair-weight reference; GenIce |
| 22 | Neighbor-distance parameter D | Implemented | `op.neighbor_distance_parameter(f_alpha=1, f_beta=1, f_gamma=1)`; powers or callables, optional symmetrization/averaging | Independent neighbor-triangle sums; GenIce serial/parallel checks |
| 23 | Angular Fourier series parameter | Implemented | `op.angular_fourier(a=1, radial=0)`; minimum-radius or product radial weights, optional averaging | Orthogonal-angle and weighted references; GenIce serial/parallel checks |
| 24 | CHILL classifier | Implemented | `op.chill()`; original narrow eclipsed band and neighbor-dependent interface rule | Ih/Ic labels; interface-rule tests; all GenIce serial/parallel checks |
| 25 | CHILL+ classifier | Implemented | `op.chill_plus()`; cutoff q3m correlations and ice/hydrate/interface labels | Exact OVITO agreement on all 11 GenIce phases and three perturbed surface configurations |

## Other methods discussed in the manuscript

Non-neural methods outside the initial 25-calculation scope are also
implemented. The optional extras keep descriptor/template dependencies out of
ordinary imports.

| Method | Kind | Status | API / evidence |
|---|---|---|---|
| Local-neighborhood A1/A2/A3 and local averages | Geometric OP family | Implemented | `op.neighborhood_parameters(m=1, averaged=False)`; independently checked sums, source-verified MALIO A/P/N convention resolves duplicated equations |
| Bond-angle B and local average | Geometric OP family | Implemented | `op.bond_angle_parameter(n1=1, n2=1, phi=0)`; mean cos(n2 theta + phi)^n1, checked against direct angle sums |
| Angle-histogram / Fourier-amplitude H | Descriptor family | Implemented | `op.angle_histogram_parameter(frequencies=[1,2,3], bins=180)`; explicit normalization; `power=True` returns MALIO-style squared amplitudes; NumPy DFT reference |
| Russo-Tanaka zeta | Geometric OP | Implemented | `op.zeta(cutoff=3.5, angle=30)`; first non-H-bonded minus last H-bonded O--O distance, H-bonds counted either way; brute-force reference on all GenIce phases |
| Ring statistics | Topological OP family | Implemented | `op.ring_statistics(max_size=10, network="hbond")`; smallest ring through every bond angle (King/Guttman), periodic wrap-around paths rejected; exact square-lattice and wrapped-chain references, known ring sets for Ih, Ic, VI, VII, VIII and XII |
| Quadrupole-related Q_T | Molecular quantity | Implemented; inputs required | `op.quadrupole_parameter(dipole, negative_site_distance=0)`; three-charge model geometry with explicit dipole/site offset; analytical test |
| Polyhedral Template Matching (PTM) | Structural classifier | Implemented; optional backend | `op.template_matching()`; install `[ptm]`; diamond-template labels verified on Ih/Ic |
| SOAP | Vector descriptor | Implemented; optional backend | `op.soap(...)`; install `[descriptors]`; serial/process, rotation, center ordering and triclinic GenIce translation tests |
| MBTR / LMBTR | Vector descriptors | Implemented; optional backend | `op.atomic_descriptor("mbtr"/"lmbtr", symbols=..., **settings)`; explicit geometry/grid/weighting; backend integration tests |
| Atom-centered symmetry functions (ACSF) | Vector descriptor | Implemented; optional backend | `op.atomic_descriptor("acsf", symbols=..., r_cut=..., g2_params=...)`; center ordering and serial/parallel tests |
| Graph-based descriptors | Descriptor family | Implemented | `op.graph_descriptor()`; directed periodic edges, topology-index mapping and invariant degree/distance features; analytical fixture |
| Coulomb, sine and Ewald-sum matrices | Matrix descriptors | Implemented; optional backend | `op.atomic_descriptor("coulomb"/"sine"/"ewald", symbols=...)`; exact Coulomb/sine references and finite Ewald output |
| DeepIce | Learned classifier | Deferred by user | Trained weights, preprocessing and model-validation work reserved for later |
| GCIceNet | Learned classifier | Deferred by user | Graph-model architecture, training/inference and validation reserved for later |
| SOAP autoencoder | Learned latent descriptor | Deferred by user | Encoder training, scaling and learned weights reserved for later |

## Execution and evidence

Per-atom numerical loops use Numba kernels and a shared joblib chunk helper.
Choose `n_jobs=1` for serial execution, a positive worker count, or `-1` for all
available workers. Thread and process backends are supported. Scalar density
needs only a volume division and does not launch parallel work.

The full suite currently passes **148 tests**, including optional-backend
checks. Validation includes 11 static GenIce structures totaling 6,696 waters.
Serial/thread agreement is tested across those structures; process agreement is additionally tested on Ih and V. Not every
implemented option has an independent reference check, as indicated above.

- [API usage, definitions and limits](order_parameters.md)
- [Extended methods and optional backends](extended_parameters.md)
- [GenIce validation report](validation/genice.md)
- [Nucleation trajectory validation: 366 snapshots, 44 sources](validation/nucleation/README.md)
- [GenIce numerical results](validation/genice_summary.csv)
- [Manuscript chapter](https://github.com/dmighty007/SolOrder/blob/main/manuscript/background/background.tex)
