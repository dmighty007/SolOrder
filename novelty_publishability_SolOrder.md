# Novelty and publishability assessment

**Work reviewed:** `manuscript/review.tex` — *Local order parameters for ice polymorphs: a unified comparison and practical guide* (D. Maity, single author), with the `SolOrder` package, the `benchmarks/` pipeline and `docs/validation/`.
**Manuscript state (27 Oct 2026):** ~15,100 words, 14 figures, 4 tables, 9 equations, 309 references (through 2026), one open `\todo{DOI after deposition}`. Abstract is 383 words.
**Literature check:** OpenAlex full-text search (6 query families), forward-citation sweeps of CHILL+ (Nguyen & Molinero, 327 citers), Shi *et al.* 2022 (26 citers) and the Doi/Takahashi series, plus a targeted web check of the closest precedent.

---

## 1. Verdict in one paragraph

The central result — a planar-interface "dead layer" δ that quantitatively predicts, and inverts, the undercounting of finite ice nuclei by structural labelers — is, as far as the searches reach, **new and publishable**. The surrounding benchmark is a genuine advance in rigor over what exists (Doi 2021; Shi 2022) but is **incremental in kind** rather than new in kind. The main obstacles are not scientific: they are (i) the hybrid review/primary-research genre, which fits no journal cleanly at the current length; (ii) the fact that the headline calibration is still partly in-sample, with the pre-registered confirmation run written but not executed; and (iii) three to five missing citations that a referee in this field will notice immediately. All three are fixable in weeks, not months.

---

## 2. Novelty, claim by claim

### 2.1 The dead-layer calibration (Sec. "A planar calibration of finite-nucleus undercounting") — **novel**

The claim: measure bulk recall *f* and δ = ∫[1 − p(z)/f]dz once on a planar ice/liquid interface (3.6 Å CHILL+, 5.5 Å PTM); then N_lab ≈ f·N(1 − δ/R)³ predicts the labeled count of a finite nucleus, and the inverse recovers N from the labeled count alone (8% pooled / 13% worst single seed, CHILL+).

I found no prior work that calibrates a labeler's interfacial deficit on a planar interface and transfers it to finite nuclei. The nearest neighbors in the literature, and how they differ:

| Prior work | What it does | Why it does not pre-empt this |
|---|---|---|
| Li *et al.*, *JCTC* **21**, 1990 (2025), "Effective Nucleation Size for Ice Crystallization" | Proposes a "corrected" critical nucleus size and a generalized nucleation theory for ice | The correction is for *defects in spontaneously grown nuclei* and CNT kinetics, not for labeler-induced surface undercounting. **Not cited — this is the single most dangerous omission.** |
| Cheng, Tribello & Ceriotti, *JCP* (2017), "…from atomistic nuclei to the planar limit" | Connects planar interfacial free energy to finite nuclei via a 1/R expansion | Same methodological spirit (planar → finite, leading 1/R term), applied to thermodynamics rather than to classification. **Not cited; should be, as the structural precedent for the moment expansion in Eq. (9).** |
| *JCP* (2025), "The bcc coating of Lennard-Jones crystal nuclei vanishes with a change of local structure detection method" | Shows the apparent surface structure of a nucleus is an artifact of the detection method | Closest in message — detection method controls what you see at a nucleus surface — but offers no correction. **Not cited.** |
| Espinosa 2016; Sánchez-Burgos 2022; Gispen 2024; Blow 2021 | Establish that nucleus size and rate depend on the order parameter | Cited. They state the problem; none quantifies or inverts it. |
| Shi *et al.*, *PCCP* (2022) | Compares order parameters on the quasi-liquid layer | Cited. Qualitative disagreement at interfaces, no size correction. |

The treatment is also unusually honest about its own status: the moment expansion (Eq. 9) is shown to be an approximation, the two competing curvature effects are named, the convention dependence (1 Å of surface definition ≈ 20% in N) is quantified, and the failure below ~200 molecules is reported. That candour is an asset with referees, not a liability.

**Assessment:** publishable as a primary result on its own. It is the part of the work I would build a paper around.

### 2.2 The unified benchmark — **substantial, but incremental**

Prior art that overlaps:
- **Doi, Takahashi and co-workers (2021 series, + CGD 2023 for clathrates)** screened 493 order parameters and combinations for Ih/Ic/liquid and the III–V–L and V–VI–L triple points. Cited.
- **Shi *et al.* 2022** compared labelers on the quasi-liquid layer. Cited.
- Several ML classifier papers (DeFever 2019, Kim 2020, Moradzadeh 2023, IceCoder 2025, score-based denoising 2024) each benchmark on their own data. Cited.

What is genuinely added here, and should be foregrounded rather than buried in a 15k-word review:
1. **Thermal vs. energy-minimized** frames, showing that local averaging is the single decisive design choice (ice/liquid 0.98–1.00 averaged vs. 0.78–0.85 unaveraged). Doi's screening was largely on relaxed structures.
2. **Threshold transfer** across five independent preparations, two water models + mW, and two temperatures — with held-out scores, not oracle scores. This is the methodological point that most of the field skips.
3. **False-nuclei testing with frozen rules** — reporting the largest spurious cluster in liquid rather than only per-molecule accuracy. Operationally important and, to my knowledge, not systematically done before for water.
4. **Neighbor-convention sensitivity scan** (cutoff, averaging radius, entropy radius) with the "published threshold reused with a different convention" test.
5. **Negative results**: no local parameter separates proton order/disorder partners (Ih/XI, V/XIII not at all); no single parameter separates every high-pressure ice from its own-pressure liquid.

Reviewers will accept 1–5 as a contribution. They will *not* accept "we computed 3×10⁷ environments" as one; the manuscript already says so explicitly, which is the right call.

### 2.3 Taxonomy and ML survey — **not novel, but competent**

Standard material, well organised, correct notation, and the bibliography is current to 2026. Its value is pedagogical. Note the asymmetry a referee will flag: the abstract contrasts the benchmark with ML classifiers that "report accuracies above 99% on their own datasets", but no learned classifier is benchmarked. Either soften the contrast or add one (IceCoder and the score-based denoising classifier both have public implementations).

### 2.4 Confined-water transfer test (supplement) — **weakest element**

Largely negative (the parameters respond to missing neighbors, not to distinct structure), on systems — host–guest OAMe–G2, T4 lysozyme — unrelated to ice, with trajectories that are *not* in the archive. It invites "out of scope" from referees and dilutes the paper's identity. I would cut it to two sentences in the Outlook, or publish it separately.

### 2.5 SolOrder as software — **unpublished, needs positioning**

No prior publication of SolOrder appears in OpenAlex. The review would become its citable reference. But pyscal, OVITO and freud already provide Steinhardt, CHILL+/PTM and Voronoi measures; a referee will ask what SolOrder adds. The defensible answer is in the manuscript already but unstated as a claim: *one* implementation covering all 22 families with a fixed settings registry, so that cross-parameter comparison is not confounded by code. Say that explicitly in the Software subsection.

---

## 3. Publishability

### 3.1 What is already strong

- **Verification before comparison**: q_l against the Legendre addition formula to 2×10⁻¹⁵, CHILL+ against OVITO (exact on all GenIce phases), mW against an independent energy calculation to 3×10⁻¹⁵, serial vs. parallel on 366 frames. Few benchmark papers do this.
- **Inferential discipline**: independent preparations (not molecules, not frames) treated as the unit of evidence; block bootstrap; explicit statement that oracle scores are not predictive accuracies.
- **Liquid-reference diagnostics**: the 220 K TIP4P/Ice liquid is explicitly shown to be slow but non-aging, and the 3.0 GPa liquid flagged as nearly arrested. Referees look for exactly this and usually do not find it.
- **Negative controls**: a bug in an existing F4 implementation found and fixed; ring-counting through periodic images flagged.
- **Reproducibility**: checksummed archive, fixed splits, a few-second quickstart that regenerates a published number, figures generated from committed CSVs.
- **Figures**: F10 is publication-quality as it stands — four panels, consistent colour coding, the model curve overlaid on individual systems and means.

### 3.2 Blockers and risks, in priority order

**(1) The confirmation run has not been executed.** `benchmarks/confirmation_protocol.md` is a proper pre-registration — frozen code hash, frozen calibration file, pass criterion stated in advance, run G seed sizes (300, 600) unused before, plus prism facets and a 500 ps time series. It was written today; no run-G results exist in `benchmarks/results/`. As the protocol itself states, the dead-layer definition and the 200-molecule validity floor were settled *after* looking at run F, so run F validates the calibration only in part. Executing run G converts the headline from "validated on held-out seeds" to "pre-registered and confirmed", and the prism/long systems answer the two questions a nucleation referee will certainly ask (facet dependence, time dependence). **This is the highest-value remaining task and I would not submit without it.**

**(2) Genre and venue.** A 15k-word, 14-figure, 309-reference manuscript that is ~55% review and ~45% primary research fits no standard slot. Review journals (*Chem. Rev.*, *Chem. Soc. Rev.*, *WIREs CMS*) are invitation-led and resist new primary data; primary journals will question the review bulk. Two viable routes:

   - **Split (recommended).** Paper A, primary: the dead-layer calibration and nucleus-size recovery, carrying run F + run G, the labeler depth analysis, the false-nuclei section and the minimum benchmark needed to support them → *JCTC* or *J. Chem. Phys.* Paper B, review/perspective: taxonomy, ML survey, transferability and sensitivity results, the practical guide and checklist → *JCP* Perspective, *PCCP* Perspective, or *Chem. Soc. Rev.* if an invitation can be obtained. Each is then a clean, strong submission, and Paper A can cite Paper B's benchmark archive.
   - **Submit whole** as a *J. Chem. Phys.* Perspective or a *PCCP* Perspective, cutting the confined-water supplement and trimming the taxonomy to the parameters actually benchmarked. Riskier; expect a "please split" from an editor.

   A single-author review of this breadth without an invitation is harder to place than the science warrants — the perspective framing, with the calibration as the stated centrepiece, is what makes it placeable.

**(3) Missing citations.** Add and distinguish, at minimum: Li *et al.* *JCTC* 2025 (effective nucleation size); Cheng/Tribello/Ceriotti *JCP* 2017 (planar limit); the *JCP* 2025 bcc-coating/detection-method paper. Also worth adding: *JCP* **160**, 5.0193340 (2024) high-dimensional OPs and NN classifiers for amorphous ices; CEGANN (*npj Comput. Mater.* 2023); *Chem. Rev.* 2025 on solid–liquid interfacial free energy; PRX 2021 nonclassical nucleation pathways in stacking-disordered crystals. The first three are the ones that could cost a revision round.

**(4) Generality of δ.** One temperature, one model, basal facets, compact spherical seeds, observed within 20 ps of release. The manuscript says this plainly, but a nucleation referee will still ask for at least facet dependence and a second temperature before accepting the correction as usable. Run G's prism slabs cover the first; a 240 K repeat of three seed sizes would cover the second cheaply.

**(5) Statistical resolution on δ.** δ rests on nine basal slabs from three preparations. Quote an uncertainty on δ itself and propagate it into the error bars of Eq. (8) — currently the errors reported are calibration-vs-truth discrepancies, not confidence intervals.

**(6) Housekeeping.** Abstract at 383 words exceeds most journals' limits (JCP 250, JCTC ~250) and currently carries five distinct results; lead with the calibration and move the rest to the introduction. Resolve the `\todo{DOI after deposition}` by depositing the archive (Zenodo) before submission. Figures are included via `../benchmarks/results/figures/` — repackage as a self-contained submission tree. The repository carries `new.tex`, `new_v1.tex` and `review.tex` as parallel drafts of the same title; prune before anyone else looks at it.

**(7) Claims to tighten.**
   - "Every labeler misses a surface layer" — three labelers were tested, none of them learned. Say "every labeler we tested".
   - "22 parameter families, 39 outputs and three labelers" is accurate and should stay; "about 3.0×10⁷ molecular environments" should never appear without the adjacent caveat it currently has.
   - The nucleation-path application compares parameters against a CHILL+-defined core, and the manuscript already notes this is agreement with a reference labeler rather than accuracy. Keep that sentence prominent; it is the kind of self-limitation that buys credibility elsewhere.

---

## 4. Recommended sequence

1. Execute run G exactly as the protocol specifies, report the result whether it passes or fails, and add the prism-facet and 500 ps time dependence of δ.
2. Add the three missing precedent citations and write the one paragraph that distinguishes this correction from Li *et al.* 2025.
3. Decide the split. If splitting, Paper A (calibration) is close to submittable once step 1 lands.
4. Deposit the archive, fix the abstract length, repackage figures, prune duplicate drafts.
5. Cut or spin off the confined-water supplement.

**Expected outcome if steps 1–4 are done:** Paper A is a credible *JCTC* / *JCP* submission with a good chance at acceptance after one round of revision; Paper B is a credible perspective. **Submitted as-is**, the most likely outcome is a request for major revision centred on (a) the in-sample element of the calibration and (b) the review/primary hybrid.

---

### Scope of this check
Literature coverage is from OpenAlex and targeted web search; it is good for indexed journal articles and arXiv, weaker for very recent preprints and for ideas buried inside papers whose titles and abstracts do not announce them. I did not re-run the benchmark, re-derive any reported number, or verify the simulation protocols beyond reading the scripts' structure; the assessment of correctness rests on the manuscript's own internal verification record.
