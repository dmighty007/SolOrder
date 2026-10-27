# Research workflow following `manuscript/new.tex`

Reviewed 26 October 2026 against the manuscript, benchmark scripts, seed-test
construction/MD scripts, and nucleation validation report. This is the follow-on
plan for the existing [benchmark design](benchmark_design.md). It specifies
proposed work; the experiments and helper modules below have not been run or
implemented by writing this document.

## Scientific focus and scope

Build the review around one question: **which local structural information
survives thermal motion, interfaces, and changes in the analysis convention?**
Compare three distinct uses of an OP: describing a distribution, assigning a
local label, and driving a simulation as a collective variable (CV). Success at
one task does not establish success at the others.

The main paper should contain the taxonomy, reproducible scalar/labeler
benchmark, thermal-versus-quenched comparison, neighbor sensitivity, cost, and
one carefully qualified nucleation application. Start with Ih/Ic/Isd and liquid,
then extend thermal coverage to selected high-pressure phases. Retain the
18-phase ideal comparison as a separate structural reference.

DeepIce, GCIceNet, autoencoders and learned CVs remain literature discussion
and future work, following the decision to defer complicated ML. A transparent
two-feature logistic baseline is sufficient for the first complementary-OP
study. SOAP/ACSF/LMBTR comparisons can follow after the scalar baseline is sound.
Amorphous ice is optional: either supply a validated preparation protocol and
data, or remove the benchmark promise from the abstract and discuss it in the
review/outlook. Empty clathrate frameworks must be identified as empty frameworks.

## What the current evidence supports

| Evidence | Present status | Limit on interpretation |
|---|---|---|
| D1: ideal GenIce crystals | 18 phases, five generator seeds, 89,100 environments; 39 scalar outputs and three labelers in the registry | Generator geometry, coordinate precision and 0.02 Å atomwise jitter affect rankings; ordered structures repeat across generator seeds |
| D4: thermal/quenched seed test | 12 systems, ten thermal and ten matched quenched frames each; 834,040 environments | One short run per system; frames 1 ps apart; construction labels require a survival audit |
| Current F1/F2/F7 | Separability, pooled redundancy and seed-label scores already generated | Thresholds are optimized on the evaluated data; frame bootstrap does not account for temporal blocks or independent replicas |
| Nucleation validation | 366 sampled frames from 44 sources; numerical and serial/parallel checks documented | CHILL+ core labels are not independent truth; recovered ordering and approximate NPT boxes restrict some sources to descriptive use |
| D2/D3 | Run E: six high-pressure ices and two compressed liquids at 220 K, two preparations each; run D: TIP4P/2005 at 220 and 200 K with a TIP4P/Ice protocol control, and oxygen-only mW at 220 K | One temperature and pressure per phase; mW at one state point; no three-site model |
| D6 | Amorphous pilot (`benchmarks/md_amorphous.py`, TIP4P/Ice, one preparation): Ih compressed at 80 K amorphizes at 1.3-1.55 GPa; recovered HDA 1.19 g/cm3 at 1 bar (CHILL+ ice 0%, CN 6.3, d5 3.25 A); liquid hyperquenched at 10 K/ns gives an LDA-like glass of 0.950 g/cm3 at 80 K (CN 4.1, d5 3.8 A, Q 0.88, CHILL+ ice 1.4%); heating HDA at 1 bar, 10 K/ns, relaxes gradually to 0.985 g/cm3 at 200 K: no distinct HDA→LDA transition under this 10 K/ns protocol | One preparation, fast rates; not yet benchmark data or in the paper |

The seed MD script uses 20 ps of preparation, with seed oxygens restrained
during preparation, followed by ten frames over the first 10 ps after release.
The velocity initialization seed is fixed in that script. This is a useful
short-time construction test, not a converged ensemble of independent nucleation
or equilibrium trajectories. New replicas must vary initialization deliberately.

## Execution order and decision gates

| Stage | Work | Required output | Gate before advancing |
|---|---|---|---|
| W0 | Freeze the pilot and audit claims/settings | Dataset manifest, claim-to-result table, applicability table | Every plotted result can be traced to coordinates, box, labels, code and settings |
| W1 | Correct statistical evaluation | Grouped splits, paired block uncertainty, separate descriptive and held-out scores | No molecule-level leakage; inadequate replication is explicitly marked |
| W2 | Audit seed truth and explain labeler failures | Survival/strain audit, interior/interface breakdown, bond-correlation diagnostics | Construction membership and current structural state are distinguished |
| W3 | Neighbor, averaging and finite-size sensitivity | F5 and uncertainty on ranking stability | Default results preserved; all tuning confined to development data |
| W4 | Generate converged thermal references | Extended F1/F4, thermal–quenched pairs, phase-retention reports | Verified state points, cell treatment, equilibration and independent replicas |
| W5 | Test information complementarity | Revised F2 and held-out two-OP maps F3 | Improvement over the best single OP persists across held-out runs |
| W6 | Measure cost and assess CV suitability separately | F6 and an evidence-based CV table T3 | Performance measurements reproducible; force-ready claims require gradients and boundary tests |
| W7 | Apply frozen methods and finish the guide | Nucleation comparison F8; recommendations T2; reproducible manuscript | Application trajectories untouched by calibration; claims match available datasets |

W1 and W2 use existing data first. W3 can begin on the audited pilot, then be
confirmed on W4. Prepare W4 inputs after W0–W2; expensive production waits for
the simulation pilot to pass. W5 needs the replicated W4 data for confirmatory
results. Cost measurements in W6 can proceed independently after settings are
frozen. The differentiable-CV development branch is optional and need not delay
publication of a qualified CV suitability assessment.

## W0–W1: make the benchmark statistically reviewable

1. **Freeze provenance.** Save code commit or source-content hashes, dependency
   versions, registry/configuration hashes, all input hashes, units, selections,
   generator/MD/noise seeds, frame times and per-frame periodic cells. Include
   hashes in cache keys so changing settings cannot silently reuse stale arrays.
   Assert that evaluated coordinates and cells match the requested snapshot.
2. **Separate outputs from methods.** Report 39 scalar columns rather than 39
   independent methods: one registry entry can emit several correlated columns.
   Record hydrogen requirements, neighbor requirements, missing-value rates and
   applicability for every condition. Define common eligibility masks for direct
   comparisons, and also report coverage on the full dataset. Do not let silent
   finite-value filtering improve a method's apparent score.
3. **Keep KS as a descriptive metric.** For two empirical distributions,
   `(1 + D_KS)/2` equals the best balanced accuracy achievable by a single
   threshold with either direction on those same observations. Call this an
   *oracle single-threshold score*. It is not held-out predictive accuracy.
   Report `D_KS` and this score separately from test-set balanced accuracy.
4. **Define experimental groups before training.** Hold out complete independent
   replicas or trajectories, with phase coverage in every split. For a single
   sufficiently long trajectory, use separated contiguous time blocks and state
   that the result measures within-run temporal transfer. Never random-split
   molecules. All threshold direction, scaling, feature selection and tuning
   occur inside training/validation groups; the outer test groups stay frozen.
   This follows the rationale for
   [grouped cross-validation](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation-iterators-for-grouped-data).
5. **Estimate uncertainty at the sampling level.** Measure autocorrelation of
   frame-level summaries, including density, energy, structural measures and
   seed size where applicable. Use blocks long enough to cover relevant slow
   correlations and check sensitivity to block length. Resample independent
   runs, then time blocks within runs. Preserve shared block IDs for paired
   thermal/quenched comparisons and hexagonal/cubic populations from one frame.
   Report runs, blocks and effective sample size, not only molecule counts.
6. **Do not manufacture precision.** Ten closely spaced frames from one run
   cannot establish between-run uncertainty. Keep present numbers as pilot
   estimates where replication is absent. Use a fixed bootstrap seed and at
   least 1,000 resamples for final intervals, checking interval stability; more
   resamples cannot compensate for inadequate independent data.

**Status (26 October 2026, partial).** Oracle scores are now labeled as
such, and held-out scores are computed against a second existing run (run B:
`train_*.npz` from `md_train.py`, 60 frames 5 ps apart, bulk Ih/Ic/Isd/liquid
at 220 K and Ih/Ic at 250 K). Results are in
`benchmarks/results/separability_heldout.csv`: held-out thresholds lose at
most 0.011, so threshold optimism is negligible at these sample sizes; the
choice of liquid reference changes distance/density-based scores by up to
0.05. Update: five independent preparations per condition (run C1–C5,
`benchmarks/md_replicas.py`; new proton configurations, new Isd stacking,
re-melted liquids) give a median between-preparation SD of 0.001 and
leave-one-preparation-out losses ≤ 0.001 wherever the oracle score is ≥ 0.8;
separations below ~0.65 can be preparation-specific
(`separability_grouped.csv`). Second update: a second water model (run D,
TIP4P/2005 at the same temperature and the same undercooling, with a
TIP4P/Ice control on the same NPT protocol) and thermal high-pressure phases
(run E, II, III, V, VI, VII, VIII and liquid at 0.3, 0.5, 1.0 and 3.0 GPa) are in
`separability_model.csv` and `separability_pressure.csv` (W4, below). Still
open: the provenance items in step 1.

Proposed implementation: extend `benchmarks/metrics.py` and add small helpers
for dataset/split metadata and hierarchical resampling. Add meaningful synthetic
checks for group separation, paired resampling and train-only threshold fitting
when implementing these helpers. Preserve pilot CSVs under a clearly versioned
result directory before regenerating final figures.

## W2: resolve the small-seed result before recommending labelers

**Status (27 October 2026).** Run F (`benchmarks/md_interface.py`,
`interface_depth.py`, F10): replicated seeds of 50–1600 molecules (three
preparations, Ih/Ic/Isd) and planar slabs. Recall against depth depends on
curvature; the planar dead layer predicts labeled counts of seeds of 200–1600
molecules via N(1 - delta/R)^3 within 9% (CHILL+) and 18% (PTM). False nuclei
in held-out liquids with frozen rules are in `false_nuclei.py` (F9).

The reported CHILL+ correct-polytype recall drops from nearly complete bulk
labeling to about 11% for the 50-molecule Ih seed interior. Treat this as a
scientific question to resolve, not yet as an intrinsic size limit of a method.

- Keep construction ancestry, construction stacking and current structural
  state as separate fields. Audit reference-neighbor retention, displacements,
  local strain and stacking changes using the known constructed lattice.
  Compare thermal snapshots with their matched minimized snapshots. Flag
  ambiguous or reconstructed environments rather than silently calling them
  classification errors or removing them using the OP being evaluated.
- The current interior flag means all four original lattice neighbors are
  inside the seed. That may still truncate a second shell or an averaging
  neighborhood. Report distance to the surface and complete-neighborhood
  coverage as a function of each method's spatial support.
- Report ice detection and correct Ih/Ic assignment separately. Include
  precision, recall, liquid false positives, wrong-polytype fraction, abstention
  or other labels, and detected versus constructed cluster size. Show core,
  interface and far-liquid populations separately.
- Plot CHILL/CHILL+ bond-correlation distributions and eclipsed/staggered bond
  counts. Test whether the eclipsed window explains the bulk Ih difference.
  Audit PTM RMSD distributions versus depth and strain. Any threshold tuning
  uses development replicas only; retain literature-default results alongside it.

**Status (26 October 2026).** Done in `benchmarks/seed_audit.py`, figure
F7b. Survival: 81–90% of molecules in 100/200-molecule seeds are intact
(within 1 Å of the constructed site after drift removal, all constructed bonds
< 3.5 Å); the 50-molecule Ih seed is the most deformed (57%). Recall is a
function of depth below the constructed surface for all sizes and polytypes:
CHILL+ ≤15% within 3 Å, ~70% at 4–5 Å, >90% beyond 6 Å; PTM needs ~2 Å more.
CHILL's bulk-Ih deficit is confirmed: the eclipsed correlation falls in its
window for 70.8% of molecules (CHILL+ window: 100%). Not yet done: local strain
and stacking change per molecule, and a development-replica threshold scan.

Gate: the F7 caption states exactly what its reference labels mean, and the
small-seed conclusion separates interface/support effects, real deformation
and classification errors. No evaluated labeler supplies its own reference truth.

## W3: robustness of the structural information

**Status (26 October 2026, partial).** `benchmarks/sensitivity.py`, figure F5:
first-shell cutoff 3.2–3.8 A for every cutoff-based column and CHILL/CHILL+,
LSI shell 3.2–3.8 A, averaging radius 3.2–5.0 A (q_lm and Q at 3.5 A) and
entropy radius 4–6 A, one factor at a time on C1–C3 at 220 K, with
cross-convention held-out scores. Leading parameters move by ≤ 0.02 with the
first-shell cutoff and their thresholds transfer within 0.025; unaveraged
first-shell parameters do not (CN 0.53–0.81). The averaging and entropy radii
matter more than the cutoff (q3 averaged Ih/Ic 1.000 → 0.65 at 4.5 A; the
perfect quenched Ih/Ic entropy separation exists only at r_m = 5 A). Not yet
done: coordinate resolution, cell size, ring definitions, other
temperatures and the high-pressure phases.

Use canonical definitions as the baseline and assign new registry IDs to
experimental variants. The following are proposed development grids, to be
frozen before confirmatory evaluation:

| Factor | Experiment | Important distinction |
|---|---|---|
| First-shell rule | Applicable radial cutoffs 3.2–3.8 Å in 0.1 Å steps; compare supported nearest-neighbor variants | Compare definitions explicitly; a different neighbor rule changes the OP |
| Averaging support | Unaveraged and independently varied averaging radii, e.g. 3.5/4.5/5/6 Å | Separate denoising gain from loss of interface resolution; some APIs currently couple these radii |
| Search completeness | Increase search radius and record missing kth neighbors/CSP12 coverage | LSI's physical shell radius (3.7 Å by default) differs from its search radius, which must exceed the shell |
| Smooth coordination | Vary tail truncation independently of `r0`; inspect boundary crossings | A decaying switch with a nonzero value at the hard cutoff is still discontinuous there |
| Ring statistics | Maximum ring size 6/8/10, cell replication, periodic winding audit | State the exact ring definition and denominator; distinguish contractible rings from periodic loops |
| Coordinate resolution | Zero jitter, high-precision coordinates and 0.01/0.02/0.05 Å perturbations | Independent atom jitter changes O–H geometry; compare rigid-molecule perturbations for hydrogen-dependent methods |
| Cell size | At least two adequate supercell sizes for shortlisted methods | For triclinic cells use minimum cell height, not just edge length; all search radii must fit, followed by convergence checks |

For proton-order partners, add a controlled comparison with identical oxygen
coordinates and cell but different valid hydrogen arrangements. Separately test
the native relaxed lattices. This distinguishes sensitivity to proton ordering
from sensitivity to density/strain. Failure of the tested hydrogen-aware local
summaries does not prove that all proton-order observables fail; polarization
and longer-range orientational correlations remain possible future additions.

Gate: report ranking changes and intervals across conventions, including cases
where a recommendation reverses. Do not select one convenient cutoff on the
final test set.

## W4: thermal references with a controlled simulation budget

**Status (26 October 2026).** First, second and transfer rows done for one
state point each (`benchmarks/md_conditions.py`, runs D and E): Ih/Ic/Isd/
liquid replicated five times (run C); II and III at 0.3 GPa, V at 0.5, VI at
1.0, VII and VIII at 3.0 GPa, liquid at 0.3, 0.5, 1.0 and 3.0 GPa, 220 K, two GenIce
preparations each, NPT then NVT, all phases retained; TIP4P/2005 Ih/Ic/Isd/
liquid at 220 K and Ih/Ic/liquid at 200 K on the C1–C3 structures. Thermal
frames transfer across models within 0.033; quenched density-like parameters
do not (down to 0.47). Every crystal pair stays separable by some parameter;
compressed liquid against II, III and V is the hard case (0.85–0.97). mW
(`benchmarks/md_mw.py`, OpenMM Stillinger-Weber, 220 K, C1–C3 oxygens) keeps
the leading parameters and their thresholds within 0.01, but first-shell
distance parameters fail to transfer even in thermal frames (smooth CN
0.90 → 0.67), and quenched mW Ih/Ic separate perfectly through lattice
density alone. The third-priority phases, other state points, a three-site
model and amorphous ice remain open.

| Priority | Conditions | Purpose |
|---|---|---|
| First | Replicated Ih, Ic, Isd and bulk liquid at the existing 220 K reference condition; selected additional temperatures | Replace short pilot sampling; distinguish bulk denoising from interface performance |
| Second | II, III, V and VI plus suitable liquid references at verified model-specific state points | Test high-pressure separability beyond ideal generator geometries |
| Third | Remaining feasible crystalline phases and selected empty frameworks | Broaden coverage while explicitly documenting metastability or transformation |
| Transfer | TIP4P/2005 subset and oxygen-only mW subset | Quantify model dependence after the first-model protocol is frozen |
| Optional | Independently prepared LDA/HDA | Include only with preparation/history and structural validation |

Start with three independent preparations per condition and expand where
between-run uncertainty remains large. Treat this as a starting design, not a
guarantee of adequate sampling. Record proton configurations and stacking
realizations; independent velocities alone do not sample proton configurations.

Before production, verify each state point against the phase diagram for the
actual force field. The candidate pressures/temperatures in the original design
are placeholders, not validated production settings. Distinguish thermodynamically
stable phases from kinetically retained metastable references. Check density,
enthalpy, cell strain, diffraction/RDF and phase-appropriate lattice/network
retention. CHILL+ is not a generic high-pressure-ice stability test; selecting
frames by the OP being ranked would also introduce selection bias.

For cell relaxation, choose and validate the ensemble explicitly. OpenMM's
[anisotropic barostat](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.MonteCarloAnisotropicBarostat.html)
scales individual axes; its
[flexible barostat](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.MonteCarloFlexibleBarostat.html)
allows independent lengths and angles. For nonorthogonal crystals requiring
shear relaxation, assess the flexible option available in the installed version
or justify a constrained-cell protocol with strain controls.

Use the old 2 ns equilibration/10 ns production proposal only as a budget
starting point. Set final lengths and analysis spacing from equilibration,
autocorrelation and between-replica convergence. Save native per-frame cells
and times. Maintain the same coordinates/frame IDs for thermal–quenched pairs;
record minimization tolerance, constraints and whether the cell is fixed.

Compare phases and liquid at matched T/P where possible, explicitly flagging
metastable liquids. Across models state whether matching uses absolute T/P or
reduced undercooling, and use a phase-appropriate reference. mW differs in
molecular resolution as well as parametrization; compare only compatible
oxygen-based observables. Measure pilot ns/day and memory on available hardware
before quoting a completion date for the thermal campaign.

## W5: complementary OPs and robust recommendations

Recompute redundancy both within each phase/condition and on a balanced pooled
sample. Examine whether pooled correlations simply encode ice versus liquid or
density. Check dependence on temperature, pressure and run. A dendrogram cut
does not establish universal equivalence between its members.

For the 39 current scalar columns there are **741 unordered pairs**. Compare
the best single OP with two-OP logistic baselines for ice/liquid, Ih/Ic and the
high-pressure subset. Fit scaling and any model settings inside grouped
validation, select pairs there, and evaluate only on held-out groups. Keep
linear and any explicitly chosen nonlinear baselines separate so feature
information is not confused with classifier capacity.

Report improvement, uncertainty and frequency of selecting each pair across
outer splits. A pair earns a recommendation only if its gain survives held-out
runs and reasonable neighbor conventions. High-dimensional descriptors and
complex trained models need not be added to make this comparison useful.

## W6: cost measurements and a distinct CV assessment

**Status (26 October 2026, cost only).** `benchmarks/cost.py`, figure F6:
every registry entry from a fresh snapshot, 1,920–65,536 molecules, 1 and 8
threads, i7-12700F. Linear scaling; whole registry 141 µs per molecule
serial (124 after the A1-A3 loop moved to a Numba kernel). Eight threads give
1.6x for the whole registry at 65,536 molecules and are slower below ~8,000. Cold import/JIT, I/O, process backends and memory are
not measured. The CV assessment has not started.

For F6, begin at roughly 10^3 and 10^4 waters; attempt 10^5 only after measuring
cost and memory, especially for topology/Voronoi methods. Report cold import/JIT,
warm calculation, neighbor construction, I/O and full-pipeline time separately.
Measure both a standalone call and batch evaluation with shared neighbors.
Use repeated timings, fixed hardware and stated dependencies. Compare serial,
threads and processes with controlled worker/backend thread counts; include
total worker memory and avoid nested joblib/Numba/backend oversubscription.
Validate numerical agreement before interpreting speedup. Cross-package checks
must match neighbor rules, normalization, precision and thresholds.

Replace the registry's single `smooth` interpretation with separate evidence:
coordinate continuity, differentiability at neighbor exchanges, cutoff handling,
available gradients, and validation for biasing. Current `cn_smooth` hard-truncates
an infinite-tail switch at the search cutoff; entropy and local averaging also
use hard neighbor membership. Gaussian broadening alone does not make the
implemented observable globally differentiable. Conversely, a kth-neighbor
distance can remain continuous at a rank exchange while its derivative changes.

Assess coordinate sweeps through cutoff crossings, neighbor ties, periodic
wrapping and angular thresholds. For optional force-ready variants, introduce
explicitly versioned compact-support weights with suitable derivative continuity,
weighted averaging and verified gradients. Compare analytic/autodiff gradients
with finite differences before any force adapter or biased dynamics. PLUMED's
[switching-function documentation](https://www.plumed.org/doc-v2.10/user-doc/html/switchingfunction.html)
provides a useful implementation reference for finite-range switching, but each
chosen variant still needs its own derivative checks. Numba acceleration does
not supply coordinate gradients.

Good phase separation or correlation with a CHILL+ cluster does not establish
reaction-coordinate quality. A committor study is a separate follow-up: define
independent liquid/crystal basins, use unbiased shooting with valid native boxes
and randomized velocities, and quantify shooting uncertainty. Existing biased
paths alone do not supply an unbiased committor or nucleation rate.

## W7: nucleation application and manuscript completion

Freeze calibration before using the application paths. Prefer readable native
fixed-box E/F/I/M sources from the
[nucleation audit](validation/nucleation/README.md), retaining untouched paths
for application. Keep recovered G/J/K ordering and approximate-box O results as
explicitly descriptive controls; do not use them for confirmatory temporal or
NPT claims until exact metadata are recovered.

Compare CHILL+, PTM and calibrated scalar/pair methods with the same declared
cluster-connectivity convention. Show strict core and interface-inclusive sizes
separately, disagreements near the interface, and per-path variation. Describe
CHILL+ agreement as agreement, not independent accuracy. Construction-based seed
tests support labeling claims only after W2. Use full native sampling/times
before studying cluster persistence or growth rates; twelve snapshots per path
are insufficient for those temporal analyses. Biased/selected paths remain
structural applications rather than rate estimates.

Update the practical guide by task, listing required inputs, spatial support,
default/tuned settings, measured cost, held-out evidence and failure conditions.
Tie every recommendation to conditions actually tested. Replace abstract promises
with completed coverage, qualify pilot-only results, and update the statement
that a common implementation removes all cutoff/model confounds: standardized
declared settings make comparison reproducible, while sensitivity still has to
be measured.

Complete the taxonomy and phase summary in parallel with analysis. Update the
2023–2026 literature from primary sources, separating experimental phase facts,
model-specific behavior and published ML results. Finish with an archived
environment, dataset manifest, figure-generation commands and data/code DOI.

## First implementation sequence

| Order | Concrete change | Acceptance check |
|---|---|---|
| 1 | Version pilot results; add replica/time/block/pair IDs and richer registry applicability/CV metadata | Every observation has a declared sampling unit; cache provenance invalidates on relevant changes |
| 2 | Separate KS oracle scores from held-out threshold scores; implement grouped splits and paired hierarchical bootstrap | Synthetic leakage/pairing checks pass; current short data are labeled preliminary |
| 3 | Add seed-survival and neighborhood-support diagnostics | Re-express F7 with audited reference labels and separate ice/polytype errors |
| 4 | Add an independent first-shell/averaging/search-radius sweep on a bounded pilot subset | F5 can be reproduced without changing canonical registry defaults |
| 5 | Prepare one replicated thermal pilot and benchmark throughput | State point, phase retention, cell treatment, autocorrelation and resource estimate documented before scaling |

These are the next work packages. MD wall time and the final submission schedule
depend on the thermal pilot and available GPUs; the current short-run timings
should not be extrapolated into a campaign ETA.
