# Order-parameter benchmark

Pipeline for Chapter 5 of the review (`manuscript/new.tex`). The design is in
[docs/benchmark_design.md](../docs/benchmark_design.md). This directory
currently implements **milestones 1–2**: existing data plus generated ideal
crystals, no new simulations.

The follow-on [research workflow](../docs/research_workflow.md) sets the next
priorities: statistical evaluation, seed-label audit, sensitivity, replicated
thermal references, and qualified nucleation/CV comparisons.

Datasets:

- **Ideal crystals:** 18 GenIce2 phases (Ih, Ic, II, III, IV, V, VI, VII,
  VIII, IX, XI, XII, XIII, XIV, XVII, and empty sI, sII and sH frameworks),
  each from 5 GenIce seeds and jittered by 0.02 A. Ice XVI is left out
  because GenIce gives it the same oxygen positions as empty sII. Generate
  with `scripts/generate_genice_fixtures.py --phases ... --output
  cache/genice/seed<N>` (see `compute.py`).
- **Thermal and quenched TIP4P/Ice at 220 K:** the seed test in
  `ice_nucleation_pg` (bulk Ih, Ic, stacking-disordered ice, and seeds of
  50–200 molecules in liquid, with per-molecule ground truth). Run A.
- **Replica runs (held-out data):** separate TIP4P/Ice NVT runs
  (`md_train.py` in `ice_nucleation_pg`) of bulk Ih, Ic, Isd and bulk liquid
  at 220 K, and Ih and Ic at 250 K; 60 thermal and 60 quenched frames 5 ps
  apart. Run B. Same starting structures (and proton configurations) as run A.
- **Independent preparations (C1–C5):** `md_replicas.py` (OpenMM + CUDA, e.g.
  the `gencompas` conda env; about 25 min for all 30 runs on an RTX 3080).
  New proton configurations (ice-rule loop moves, net dipole < 2%; the run A/B
  crystals have about 18%), new stacking sequences for Isd, liquid re-melted
  at 300 K and re-equilibrated 1 ns at 220 K. Same frames and protocol as run B.
- **Second water model (run D):** `md_conditions.py`. TIP4P/2005 on the C1–C3
  structures: Ih, Ic, Isd and liquid at 220 K (same temperature) and Ih, Ic,
  liquid at 200 K (same undercooling; melting points about 252 and 272 K), plus
  a TIP4P/Ice control with the same protocol. NPT at 1 bar (anisotropic for
  crystals; 200 ps crystals, 1–2 ns liquids), then NVT in the final box with
  the run C frame schedule. Separates the protocol effect (C → TIP4P/Ice NPT)
  from the model effect (TIP4P/Ice NPT → TIP4P/2005).
- **Oxygen-only mW (run D):** `md_mw.py`. mW as an OpenMM Stillinger–Weber
  potential (CustomNonbondedForce + CustomManyParticleForce; energy equal to a
  direct numpy sum to 3e-15, CPU/CUDA/Reference agree), on the C1–C3 oxygen
  positions and the seed-test liquid at 220 K (mW melts at 274.6 K, so this
  is also the TIP4P/Ice undercooling). Same protocol as above, 5 fs steps.
  Run on solaris node1 (3 x RTX 2080 Ti, `cvs` conda env, OpenMM 8.2 CUDA
  11.8), about 4 min for all 12 runs. Hydrogen-dependent columns are NaN.
- **Thermal high-pressure phases (run E):** `md_conditions.py`. TIP4P/Ice at
  220 K: ice II and III at 0.3 GPa, V at 0.5, VI at 1.0, VII and VIII at 3.0
  GPa (GenIce seeds 1 and 2, every periodic height ≥ 24 A, 768–1500
  molecules) and liquid at 0.3, 0.5, 1.0 and 3.0 GPa (4096 molecules; each crystal has a pressure-matched liquid), two preparations
  each. NPT (flexible barostat for triclinic II and V), then NVT as above.
  Paired with run C1–C2 Ih, Ic and liquid for the ambient phases.

## Run

```sh
cd benchmarks
<gencompas>/bin/python md_replicas.py all   # run C trajectories -> cache/replicas/
<gencompas>/bin/python md_conditions.py all # runs D and E (about 75 min) -> cache/conditions/
<gencompas>/bin/python md_mw.py prepare     # mW inputs -> cache/conditions/mw_inputs/ (copy to the GPU host)
MW_PLATFORM=CUDA MW_DIR=. python md_mw.py run all   # on the GPU host; copy D_mW_*.npz/.json back
../.conda-env/bin/python compute.py     # about 5 min on 10 workers -> cache/*.pkl (not committed)
../.conda-env/bin/python metrics.py     # about 17 min, mostly the frame bootstrap -> results/*.csv
../.conda-env/bin/python metrics.py --steps heldout   # one step only (ideal, thermal, heldout, labels)
../.conda-env/bin/python seed_audit.py  # about 10 s -> results/seed_audit_*.csv
../.conda-env/bin/python sensitivity.py compute   # about 1 min -> cache/sensitivity.pkl
../.conda-env/bin/python sensitivity.py metrics   # about 1 min -> results/sensitivity.csv
OMP_NUM_THREADS=1 ../.conda-env/bin/python cost.py --jobs 1   # about 20 min -> results/cost_jobs1.csv
../.conda-env/bin/python cost.py --jobs 8                     # -> results/cost_jobs8.csv
../.conda-env/bin/python figures.py     # -> results/figures/*.png, *.pdf
```

`compute.py` reads the seed-test data from
`/home/dm/Dibyendu/Projects/ice_nucleation_pg/data/seedtest`. Input hashes,
the SolOrder version and the registry hash are written to
`cache/manifest.json`.

## Files

| File | Contents |
|---|---|
| `registry.py` | Every parameter, its fixed settings and its source (39 scalar columns, 3 labelers) |
| `results/columns.csv` | Column, group, whether hydrogens are needed, smoothness |
| `results/separability_ideal.csv` | KS D and oracle score (1 + D)/2 for 153 GenIce phase pairs per column |
| `results/separability_thermal.csv` | Same for Ih, Ic, liquid and Isd layers at 220 K (run A), with 95% frame-bootstrap intervals |
| `results/separability_heldout.csv` | Oracle scores on runs A and B; held-out balanced accuracy A→B (25 ps block intervals), B→A, first→second half of B, and 220→250 K |
| `results/separability_grouped.csv` | Runs C1–C5: oracle per preparation (mean, SD, min), leave-one-preparation-out and run-A-fitted held-out scores |
| `results/replica_autocorrelation.csv` | Lag 1/2/5 autocorrelation of per-frame means in run B, per system and column |
| `results/pilot_v1/` | Pilot CSVs before the held-out step (score column then called `accuracy`) |
| `results/redundancy_spearman.csv` | Spearman matrix, thermal Ih + Ic + liquid pooled |
| `results/label_scores.csv` | CHILL, CHILL+ and PTM against seed-test ground truth |
| `results/figures/F1a`, `F1b`, `F2`, `F7` | Phase-pair count matrix, thermal heat map, dendrogram, seed-test labels |
| `results/seed_audit_systems.csv` | Per seed system: intact fraction, displacement, bonds lost; labeler recall overall, among intact and deeper than 5 A |
| `results/seed_audit_depth.csv` | Per system and 1 A depth bin: intact fraction, PTM RMSD, and correct / wrong polytype / interfacial / other per labeler |
| `results/seed_audit_chill_bonds.csv` | Bulk ice: staggered counts and where the eclipsed-bond correlation falls relative to the CHILL and CHILL+ windows |
| `results/figures/F7b` | Correct polytype against depth below the seed surface |
| `results/figures/S1` | Full ideal heat map: every column x every phase pair |
| `results/figures/S2` | Held-out versus oracle scores, runs A and B |
| `results/figures/S3` | Lowest held-out score over C1–C5 versus mean oracle score |
| `results/separability_model.csv` | Run D: oracle per condition (C, TIP4P/Ice NPT, TIP4P/2005 220 K and 200 K, mW 220 K) and thresholds moved between them; `*_min` is the worst preparation |
| `results/separability_pressure.csv` | Run E: every pair of 8 crystals and 3 liquids, thermal oracle, held-out between the two preparations (mean, worse) and the ideal GenIce oracle |
| `results/figures/S4` | Run D: held-out versus oracle for the protocol and model moves |
| `results/figures/F8` | Run E: best held-out parameter per phase pair; ideal versus thermal oracle |
| `results/sensitivity.csv` | C1–C3 at 220 K: oracle, held-out (C1 → C2/C3) and cross-convention (C1 at the default → C2/C3 at the variant) scores per column, task and setting of the first-shell cutoff, LSI shell, averaging radius and entropy radius; labeler scores per cutoff |
| `results/cost_jobs1.csv`, `cost_jobs8.csv` | Wall time per registry entry and for the whole registry against system size, 1 and 8 threads |
| `results/figures/F5`, `F6` | Neighbor-convention sensitivity; computational cost |

## Findings so far

Scores below are **oracle scores**, (1 + KS D)/2: the balanced accuracy of the
best single threshold chosen on the same data. They are not predictive
accuracies, but see the held-out check, which shows they barely differ here.

- **Ideal crystals are mostly easy; three pairs are not.** 150 of 153 pairs
  have at least one parameter at an oracle score of 0.95 or better. The three
  that do not are:
  - Ih vs XI (best: smooth CN, 0.87)
  - V vs XIII (best: averaged q4, 0.83)
  - sI vs sII (best: 0.92)

  The first two are proton order/disorder partners whose oxygen frameworks
  barely differ. The other partner pairs separate only because GenIce's
  ordered cells are slightly distorted: III–IX through averaged entropy
  (the only parameter at ≥ 0.95), VII–VIII through 17 parameters and
  XII–XIV through 5. None of
  the hydrogen-aware parameters (F4, ζ, H-bond count, ring fractions)
  separates Ih from XI or V from XIII.
- **Most broadly useful on ideal crystals:** averaged entropy separates 93%
  of pairs at ≥ 0.95, then Voronoi asphericity and density (89–90%) and
  averaged Q (89%). H-bond count and A3 separate none.
- **Ih vs Ic needs neighbor-averaged or hydrogen-aware parameters.** At
  220 K the best are averaged q3 (1.00) and q4 (1.00), then F4 (0.94),
  crystalline connections (0.93), averaged q6 (0.78) and CSP (0.72). Every
  other parameter, including ζ and the ring fractions, is at 0.61 or below.
  Ring statistics cannot see stacking at all: Ih and Ic both have only
  6-rings. Averaged entropy separates Ih from Ic
  perfectly in ideal crystals but falls to 0.61 at 220 K.
- **Local averaging is what survives thermal noise.** Averaged entropy gives
  1.00 for ice vs liquid at 220 K (unaveraged: 0.85). Averaged Q gives 0.98
  (unaveraged: 0.84) and averaged q6 gives 0.99 (unaveraged: 0.78).
  Coordination number (0.62) and H-bond count (0.53) are nearly useless.
- **The 6-ring fraction is the best unaveraged parameter for ice vs liquid**
  (0.99 at 220 K), second only to averaged entropy (1.00). ζ reaches 0.84
  thermal and 1.00 quenched, behaving like d5 and LSI.
- **Quenching recovers most of the separability** that thermal motion
  removes. For example, d5 and LSI go from 0.83 to 1.00 for ice vs liquid.
- **Redundancy.** Clustering at |ρ| > 0.8 reduces the 39 columns to 22
  groups. The largest group is Q, bond angle, angular Fourier, the three
  angle-histogram amplitudes, q3, q6 and normalized w4: in this data they
  are essentially one quantity. Other groups are averaged q3/q4/q6 with
  crystalline connections; d5, LSI, entropy and ζ; averaged Q with averaged
  entropy; the 5- and 6-ring fractions; and A1 with A2.
- **Held-out thresholds lose almost nothing.** A threshold fitted on run A
  and applied to run B is within 0.008 of the run B oracle score for every
  column and pair (thermal and quenched); B→A within 0.011; first half → second
  half of B within 0.005. The best column per task is the same under both
  scores. With 10^4–10^5 molecules per population, one threshold barely
  overfits.
- **The reference populations matter more.** Ih vs Ic oracle scores agree
  between runs within 0.013. Replacing run A's far liquid with run B's bulk
  liquid changes 25 of 156 ice/liquid scores by more than 0.01, all distance-
  or density-based columns (CN, smooth CN, d5, Voronoi, bond angle, A1, A2,
  CSP, NDP), by up to 0.05.
- **Independent preparations (C1–C5) confirm it.** Median between-preparation
  SD of the oracle score is 0.001. Where the mean oracle score is ≥ 0.8,
  leave-one-preparation-out costs ≤ 0.001 and run-A thresholds cost ≤ 0.003.
  Below 0.65, separations can be preparation-specific: in quenched Isd,
  smooth CN separates hexagonal from cubic layers at 0.60 in one stacking
  sequence and 0.51–0.52 in four others (held-out 0.45). H-aware columns
  differ by ≤ 0.009 between the polarized run A crystals and C.
- **Ih vs Ic at 250 K** keeps the same order: averaged q3 (1.00), averaged q4
  (0.99), F4 (0.91), crystalline connections (0.88). Thresholds fitted at
  220 K lose at most 0.03 (averaged entropy, already weak for this pair).
- **Labelers on seeds.** CHILL+ is best on Ih and stacking-disordered
  seeds. On Ic seeds CHILL is level or slightly ahead (N=200: 0.57 vs 0.55).
  CHILL+ interior-seed scores reproduce the earlier seed-test analysis
  (e.g. Ih N=100: 0.465 vs 0.4651).
  - The original CHILL labels only 71% of thermal bulk Ih correctly
    (CHILL+: 99.7%). Cause, from the bond correlations: the eclipsed bond
    falls in CHILL's window (-0.2, -0.05) for 70.8% of Ih molecules and in
    CHILL+'s (-0.35, 0.25) for 100%.
  - PTM at OVITO's default RMSD cutoff never labels liquid as ice, but it
    misses most small seeds.
  - Liquid false-positive rates are at most 1% for all three.
- **Seed audit: small-seed recall is a depth effect** (`seed_audit.py`).
  81–90% of molecules in 100- and 200-molecule seeds are still intact
  (within 1 A of the constructed site, all constructed bonds kept); the
  50-molecule Ih seed is the most deformed (57%). CHILL+ recall rises from
  ≤15% within 3 A of the seed surface to ~70% at 4–5 A and >90% beyond 6 A,
  for every size and polytype, while ~90% of molecules at 2–4 A are intact.
  It labels 17–58% of molecules at 1–4 A as interfacial ice. A 50-molecule
  seed has only 3–4% of molecules deeper than 5 A. PTM needs ~2 A more:
  its RMSD drops below the 0.1 cutoff for half of the molecules only at
  6–7 A.

- **Second water model** (run D, `separability_model.csv`, S4). Thermal
  frames: oracle scores change by a median 0.004 (max 0.025) from TIP4P/Ice to
  TIP4P/2005 at either 220 K or 200 K, and thresholds moved across models lose
  at most 0.016–0.033 in the worst preparation wherever the oracle is ≥ 0.8
  (only smooth CN exceeds 0.02). The protocol change (constructed NVT box →
  NPT) costs ≤ 0.012. Quenched frames are different: their distributions are
  so narrow that a threshold sits on a model-specific value. Smooth CN and
  Voronoi density fall from 0.97 and 0.82 oracle to 0.47–0.65 held-out,
  MALIO A1/A2 to 0.89–0.97, and d5 (Ic/liquid, 200 K) to 0.96. Rankings and
  the best parameters (averaged entropy, q6, q3, F4) are unchanged.
- **Oxygen-only mW** (run D, S4d). The leading parameters are unchanged and
  their thresholds transfer from TIP4P/Ice within 0.01: averaged entropy
  (1.000), q6 averaged (0.998), Q averaged (0.991) and q3 averaged (0.987)
  for ice vs liquid; q3 and q4 averaged (≥ 0.998) for Ih vs Ic. Unlike
  TIP4P/2005, thermal frames now fail too for first-shell distance
  parameters, because the mW liquid packs more neighbors into the first shell
  (smooth CN median 4.76 vs 4.39): smooth CN keeps a 0.90 oracle but scores
  0.67 with the TIP4P/Ice threshold, MALIO A1 0.80 → 0.72, d5 and LSI lose
  0.04; q6 averaged Ih/Ic 0.87 → 0.81. Distance-based ice/liquid oracle
  scores rise in mW (CN 0.63 → 0.79, NDP 0.59 → 0.80, d5 0.87 → 0.92),
  Voronoi density falls (0.68 → 0.58). Quenched mW Ih and Ic become perfectly
  separable by smooth CN, Voronoi density and asphericity, entropy and NDP
  (≥ 0.998, against 0.57–0.95 in TIP4P/Ice): without proton disorder the
  quenched crystals are near-perfect lattices, and the 0.4% Ih/Ic density
  difference alone separates them. Not a stacking signal.
- **Ih/Ic density artifact under NPT.** In NPT, Ic is 0.35% denser than Ih in
  both models (0.928 vs 0.925 g/cm3 TIP4P/2005, 220 K). In quenched frames this
  alone gives Voronoi density and NDP Ih/Ic oracle scores of 0.74–0.90,
  against 0.52–0.55 in the constant run C box. That is a lattice-parameter
  signal, not a stacking signal.
- **Thermal high-pressure phases** (run E, `separability_pressure.csv`, F8).
  All phases were retained (O displacement 0.25–0.45 A over 300 ps, dipole
  cosine 0.96 including proton-ordered II and VIII; densities II 1.18, III
  1.15, V 1.26, VI 1.35, VII 1.56, VIII 1.57 g/cm3). Every crystal pair is
  still separated by some single parameter in thermal frames (best held-out
  ≥ 0.98, except VII/VIII 0.92 with q6 averaged, 1.00 when quenched with the
  angle histogram). The hard cases are high-pressure ice against dense
  liquid: V vs 0.3 GPa liquid 0.85 (Voronoi asphericity), III vs ambient
  liquid 0.94, liquids at different pressures 0.85–0.99. Per parameter,
  thermal noise is costly: of 684 crystal pair × parameter combinations with
  ideal oracle ≥ 0.95, 31% fall below 0.95 and 10% below 0.8 at 220 K; the
  quenched frames recover all but 3.5%. Losses are largest for W4/W6 (median
  0.23–0.26), NDP (0.14), q6, F4, entropy (0.06–0.08), near zero for averaged
  and network parameters.
- **Pressure-matched liquids** (run E, liquid at 0.3, 0.5, 1.0 and 3.0 GPa).
  Each crystal against liquid at its own pressure, worse held-out of the two
  preparations, thermal: VI, VII, VIII 1.00 (averaged entropy), II 0.97
  (6-ring fraction), III 0.95 (averaged entropy), V 0.84 (q4 averaged; 0.93
  quenched). Best single parameter over all matched pairs: Q averaged, 0.69.
  Liquids at nearby pressures: 0.3/0.5 GPa 0.64, 0.5/1.0 GPa 0.75 (smooth
  CN). The 3.0 GPa liquid is nearly arrested at 220 K (O displacement 0.5 A in
  300 ps, density 1.44 vs 1.56 g/cm3 for VII) and stays amorphous.
- **Quenched held-out failures across preparations.** In run E, 62 quenched
  pair × parameter scores with oracle ≥ 0.8 lose more than 0.02 between the
  two preparations (thermal: 0), mostly with ice II and VIII (d5, averaged
  entropy, smooth CN, ζ, Voronoi density). Quenched ordered crystals have
  near-delta distributions (d5 SD 0.0002–0.003 A), and preparations differ
  by the final NPT box (0.4% in density for II), so the peaks move by more
  than their width. Same mechanism as the quenched model transfer.

- **Neighbor convention** (`sensitivity.py`, F5). One factor at a time on
  C1–C3, 220 K. First-shell cutoff 3.2–3.8 A: the leading parameters move
  by at most 0.02 (q6 averaged 0.987–0.992 ice/liquid; q3 averaged
  0.983–1.000, q4 averaged 0.993–0.997, F4 0.930–0.939, crystalline
  connections 0.906–0.925 Ih/Ic), thresholds fitted at 3.5 A lose at most
  0.025, rankings hold (Spearman ≥ 0.94). Unaveraged first-shell parameters
  all rise with the cutoff (ice/liquid, thermal: CN 0.53 → 0.81, smooth CN
  0.63 → 0.87, NDP 0.57 → 0.76, bond angle 0.73 → 0.85, A1 0.76 → 0.87);
  a smooth-CN threshold from r0 = 3.5 A scores 0.50–0.69 at other r0.
  H-bond parameters and the LSI shell barely matter. Averaging radius (q_lm
  at 3.5 A): larger radii help Q averaged (→ 1.000 at 5 A), q6 averaged Ih/Ic
  (0.776 → 0.962) and q4 averaged, but break q3 averaged for Ih/Ic (1.000 →
  0.894 / 0.647 / 0.804 at 4.0 / 4.5 / 5.0 A). Entropy radius: averaged
  entropy needs r_m ≥ 5 A for ice/liquid (0.95 at 4 A); its perfect quenched
  Ih/Ic separation exists only at r_m = 5.0 A (0.53–0.77 otherwise), a
  cutoff coincidence. CHILL+ bulk Ih drops from 0.997 to 0.961 at 3.8 A
  (it needs exactly four neighbors inside the cutoff); CHILL is insensitive.
- **Cost** (`cost.py`, F6; i7-12700F). Linear in system size for every
  entry. Whole registry 141 µs per molecule on one thread (9.2 s for 65,536
  molecules); CN 1.8, PTM 2.7, CHILL/CHILL+ 5, most others 6–13 (mainly the
  5 A neighbor search), rings 30, angle histogram 21, Voronoi 20, A1–A3 19,
  entropy 18 µs per molecule. Eight threads speed up rings (3.2x),
  crystalline connections and Steinhardt (2.5x), but the whole registry is no
  faster (0.94x) and the neighborhood parameters are 4x slower: their
  per-atom Python loop (`extended._neighborhood_chunk`) holds the GIL.
  Parallelize over frames instead.

## Caveats

- Run A frames are 1 ps apart, 10 per system. The bootstrap treats them as
  independent, so intervals (median width 0.009 in D) are too narrow. Run B
  intervals use 25 ps blocks (per-frame means are nearly uncorrelated at 5 ps
  in the crystals, lag-1 up to 0.55 for liquid ring fractions) but still
  describe one run; between-run differences are the better uncertainty.
- Runs A and B start from the same structures; only C1–C5 are independent
  preparations. Runs A–C are NVT in the constructed box; runs D and E are
  NPT-equilibrated, then NVT in the final NPT box.
- The supercooled liquids keep densifying slowly: in the last half of NPT
  their volume still changes by 0.2–0.8%. Mean O displacement over the 300 ps
  of production: TIP4P/2005 1.5–1.8 A at 220 K and 0.9 A at 200 K,
  TIP4P/Ice 1.0–1.1 A at 220 K.
- The 220 K liquid is slow (~3 A per ns) and still aging: two C liquids
  drift down by 0.05–0.09 kJ/mol per molecule over 300 ps, without ice
  forming (CHILL+ ice + interfacial stays near 10%, no trend).
- `smooth` in the registry marks a smooth functional form. Smooth CN is
  truncated at the 5 A search cutoff (each neighbor there still adds about
  0.003), and entropy uses hard neighbor membership: neither is a tested CV.
- "Liquid" means molecules more than 8 A beyond the seed radius in the seeded
  220 K systems (supercooled). Ice references are separate bulk NVT crystals.
- Ideal structures are judged at 0.02 A jitter resolution, about twice the
  `.gro` coordinate precision. Proton-ordered phases are identical across
  GenIce seeds, so their 5 replicas differ only by jitter.
- Proton-order results reflect GenIce's cell geometries. Real ordered phases
  also differ in lattice parameters from their disordered partners, but by
  amounts these cells may not reproduce.
- High-pressure phases are thermal at one temperature (220 K) and one
  pressure each, two preparations, 300 ps of production. Ice VII at 220 K is
  metastable relative to VIII in experiment; it stays proton-disordered here
  only because rigid-water MD does not order protons on this time scale. No
  ice IV, IX, XII–XVII, clathrates or amorphous ice are thermalized.
- Three water models: two four-site TIP4P variants and the oxygen-only mW,
  each at one state point (mW at 220 K only). No three-site model.
- mW runs have no center-of-mass motion remover; the box drifts 2–5 A over
  300 ps (harmless for these translation-invariant parameters). The stored
  `displacement` in the 12 mW npz files includes this drift; with it removed,
  crystal oxygens move 0.43–0.46 A (md_mw.py now reports that).
- The PTM RMSD cutoff (0.1) and the crystalline-connection threshold (0.5)
  are literature defaults, not tuned.
