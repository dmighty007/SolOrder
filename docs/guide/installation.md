# Installation

SolOrder needs Python 3.10 or newer. It is built on
[MDAnalysis](https://www.mdanalysis.org/), which reads almost every simulation
format (GROMACS `.gro`/`.xtc`, LAMMPS dumps, PDB, XYZ and more).

## Install the package

```sh
python -m pip install SolOrder
```

This also installs NumPy, SciPy, pandas, Numba, MDAnalysis, joblib and freud.
The examples in this guide also use matplotlib:

```sh
python -m pip install matplotlib
```

Check that it worked:

```python
import SolOrder
print(SolOrder.__version__)
```

!!! note "The first call is slow"
    SolOrder compiles its numerical kernels with Numba the first time they run
    and caches them on disk. Expect the very first calculation to take a few
    seconds; later calls, and later Python sessions, are fast.

## Optional extras

Some methods need extra packages. Install them only if you need them:

| Extra | Command | Adds |
|---|---|---|
| `ptm` | `pip install "SolOrder[ptm]"` | polyhedral template matching (via OVITO) |
| `descriptors` | `pip install "SolOrder[descriptors]"` | SOAP, ACSF, MBTR and other DScribe descriptors |
| `examples` | `pip install "SolOrder[examples]"` | JupyterLab, matplotlib and seaborn |

## Get the example files

The guide uses small structures that ship with the source repository, in
`docs/guide/data/`:

```sh
git clone https://github.com/dmighty007/SolOrder.git
cd SolOrder/docs/guide/data
```

| File | Contents |
|---|---|
| `ice_Ih.gro` | hexagonal ice, 1920 waters, 220 K |
| `ice_Ic.gro` | cubic ice, 1920 waters, 220 K |
| `liquid.gro` | supercooled liquid water, 4096 waters, 220 K |
| `ice_Ih.xtc` | ten frames of the ice Ih run, 5 ps apart |
| `seed.gro.gz` | an ice Ih nucleus of 1600 molecules in liquid (oxygens only) |
| `seed_truth.npz` | which molecules were placed in the nucleus |

All are thermal TIP4P/Ice configurations; the massless charge site was removed
because SolOrder works with oxygen and hydrogen positions. Run the examples from
this folder, or change the file paths.

## Development install

To work on SolOrder itself:

```sh
git clone https://github.com/dmighty007/SolOrder.git
cd SolOrder
python -m pip install -e ".[dev]"
python -m pytest
```

To build this documentation locally:

```sh
python -m pip install -e ".[docs]"
mkdocs serve
```
