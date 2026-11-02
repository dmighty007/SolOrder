# Research examples

Install SolOrder and optional notebook dependencies from the repository root:

```sh
python -m pip install -e ".[examples]"
```

`notebooks/` contains the original exploratory notebooks. These are research
artifacts, not an automated test suite. Some cells reference simulation inputs
on the original author's machine; supply your own trajectories and update those
paths before running them. From a notebook, bundled sample files are available
under `../data/`. Existing outputs have been preserved under `results/`.

The water analyses expect oxygen names `OW` and hydrogen names `HW1`/`HW2`,
with periodic dimensions supplied by the MDAnalysis Universe. The bundled XYZ
files alone do not provide all of that topology and box information.
