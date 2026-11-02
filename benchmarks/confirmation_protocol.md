# Confirmation protocol for the nucleus-size calibration (run G)

Written on 2026-10-27, before any run G system had been simulated or
analyzed. Its purpose is to test the calibration on data that played no part
in choosing any definition, convention or validity range.

## Why

The dead-layer definition and the validated range were settled after the run F
seeds had been examined: the molecular-plane versus equimolar convention, the
bin-filling rule of the integral definition, the count definition that replaced
it, and the lower limit of 200 molecules. Run F therefore validates the
calibration only in part. Run G is a fresh test.

## Frozen analysis

- Code: `nucleus_size.py` (sha256 `acd3905fd8565252080a5af77d0b2d1153487f971ccb0fdf33db560d8bb5d6a1`).
- Calibration: `results/nucleus_size_calibration.csv` (sha256
  `f3f2ed941a8c6be22d0b053df0f549e8429215419c1456cc6bda8557e6eeab13`). This is f and the count
  dead layer `delta = Gamma/n` from the nine run F basal slabs only, for CHILL+ and PTM. Neither is
  changed after this point.
- Primary test: eq. 8, the inverse, applied to the correctly labeled count and to the
  largest Ih/Ic cluster of every run G seed. The quantity compared is the constructed size,
  with the intact size as a secondary comparison.
- Reported for every seed, with no exclusions: signed and absolute errors per system and the
  mean per size.
- Pass criterion, stated in advance: the run G errors lie within the range observed for
  run F seeds of 200–1600 molecules under the count definition. These are, for the
  correct-polytype count, CHILL+ at most 9.8% and PTM at most 24.7% per system; for the
  largest cluster, CHILL+ at most 18.8% and PTM at most 28.1%. A failure is reported as one.

## Run G systems (`md_interface.py`, preparation p4, new liquid)

- `G_seed_{Ih,Ic,Isd}_{300,600}_p4`: seed sizes not used before. Same construction and
  protocol as run F: restrained for 20 ps, then 10 frames 2 ps apart.
- `G_prism_{Ih,Ic,Isd}_p4`: planar slabs with primary prism {10-10} faces, standard protocol.
  These measure the facet dependence of f and delta. They are not used to recalibrate.
- `G_long_{Ih,Ic,Isd}_p4`: basal slabs, frames at 2–20 ps every 2 ps, then every 20 ps to
  500 ps. They measure the time dependence of delta, and are not used to recalibrate. Their
  first ten frames, which follow the standard protocol, also give a fourth basal-slab
  preparation, reported as a replicate.

## Analysis commands

    python interface_depth.py compute --confirmation
    python confirmation.py
