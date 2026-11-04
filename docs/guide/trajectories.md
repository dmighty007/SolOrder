# Working with trajectories

So far every example analysed a single frame. A simulation has many. SolOrder
analyses one frame at a time; you loop over the trajectory with MDAnalysis.

## Load a trajectory

Give MDAnalysis a topology (`.gro`) and a trajectory (`.xtc`):

```python
import numpy as np
import MDAnalysis as mda
from SolOrder import OrderParameter

u = mda.Universe("ice_Ih.gro", "ice_Ih.xtc")
print(len(u.trajectory), "frames")
```

```text
10 frames
```

## Loop over frames with `refresh`

Create the `OrderParameter` once, then call `refresh(frame)` to move it to a
new frame. `refresh` copies the new coordinates and clears every cached result.

```python
op = OrderParameter(u, cutoff=5.0)

for ts in u.trajectory:
    op.refresh(ts.frame)
    q = op.tetrahedral_order()
    hexagonal = np.mean(op.chill_plus().labels == 1)
    print(f"{ts.time:5.1f} ps   <q> = {q.mean():.3f}   hexagonal = {hexagonal:.3f}")
```

```text
  0.0 ps   <q> = 0.957   hexagonal = 1.000
  5.0 ps   <q> = 0.957   hexagonal = 0.999
 10.0 ps   <q> = 0.956   hexagonal = 0.999
  ...
 45.0 ps   <q> = 0.956   hexagonal = 0.999
```

!!! warning "The object is a snapshot"
    `OrderParameter` copies the coordinates when it is created or refreshed.
    Moving the Universe on its own (`u.trajectory[5]`) does **not** update
    `op`: call `op.refresh(5)`.

## Collect results in a table

A pandas DataFrame with one row per molecule and frame is convenient for
later analysis:

```python
import pandas as pd

rows = []
for ts in u.trajectory:
    op.refresh(ts.frame)
    rows.append(pd.DataFrame({
        "time_ps": ts.time,
        "atom": op.atom_indices,
        "q": op.tetrahedral_order(),
        "q6_avg": op.steinhardt(6, averaged=True, cutoff=3.5),
        "label": op.chill_plus().names,
    }))
table = pd.concat(rows, ignore_index=True)
print(table.groupby("time_ps")[["q", "q6_avg"]].mean().head(3))
```

```text
                q    q6_avg
time_ps
0.0      0.956822  0.478184
5.0      0.957463  0.479246
10.0     0.955840  0.475993
```

Save it with `table.to_csv("order.csv")` or `table.to_parquet(...)`.

## Analyse only some frames

Slice the trajectory to skip equilibration or thin out frames:

```python
for ts in u.trajectory[2::4]:      # frames 2, 6
    op.refresh(ts.frame)
    ...
```

## Analyse only some molecules

The `selection` argument restricts the analysis, for example to molecules near
a surface. Neighbors are then searched **only among the selected atoms**, so
select generously and keep the subset afterwards:

```python
op = OrderParameter(u, cutoff=5.0)
q = op.tetrahedral_order()
upper_half = op.pos[:, 2] > u.dimensions[2] / 2   # op.pos: oxygen positions (Å)
print(q[upper_half].mean())
```

## Run in parallel

Large systems run faster with several workers:

```python
op = OrderParameter(u, cutoff=5.0, n_jobs=4)    # 4 threads; n_jobs=-1 uses all cores
```

`n_jobs=1` (the default) runs serially. Threads share the snapshot, so memory
does not grow. For a few thousand molecules the serial code is already fast,
and extra workers may not help; parallelism pays off for tens of thousands of
molecules and expensive parameters.

Next: [tips and common mistakes](tips.md).
