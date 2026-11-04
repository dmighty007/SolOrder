# Choosing a parameter

SolOrder implements dozens of order parameters. Which one should you use? It
depends on the question. This page summarizes the recommendations of the
SolOrder benchmark (see the [benchmark design](../benchmark_design.md)), then
shows how to test a choice on your own data.

## Recommendations

| Question | Good choices | Avoid |
|---|---|---|
| Ice or liquid? | local entropy $\bar s$, averaged $\bar q_6$, averaged tetrahedral order $\bar q$; six-ring fraction if you want no averaging | coordination number, hydrogen-bond count |
| Ice Ih or ice Ic? | averaged $\bar q_3$, $\bar q_4$; F4 if hydrogens are available | ring statistics, unaveraged parameters |
| Polytype labels in a nucleus | CHILL+ (Ih, stacking-disordered), CHILL or CHILL+ (Ic) | CHILL for Ih; comparing nucleus sizes from different labelers |
| Finding a nucleus in liquid | $\bar s$, $\bar q_6$, CHILL+, template matching, followed by a largest-cluster step | smooth CN, $d_5$, $\bar q_4$; counting CHILL+ `interfacial_ice` as ice |
| Telling crystalline ice phases apart | $\bar s$, smooth CN, Voronoi asphericity, $\bar q$ | $\hat w_4$, $\hat w_6$, hydrogen-bond count |
| Proton-ordered or disordered? | no local parameter tested is reliable | — |
| Oxygen-only models (mW) | $\bar s$, $\bar q_6$, $\bar q$, $\bar q_3$ | distance-parameter thresholds taken from other models |

These recommendations come from TIP4P/Ice, TIP4P/2005 and mW water at 200–250 K,
plus high-pressure ices at 220 K. Treat them as a starting point and check
them at your own conditions.

## The calls

| Symbol | SolOrder call |
|---|---|
| $\bar s$ | `eFingerprint(u, cutoff=5.0).localEntropy()` |
| $\bar q_3$, $\bar q_4$, $\bar q_6$ | `op.steinhardt(l, averaged=True, cutoff=3.5)` with `l` = 3, 4, 6 |
| $\hat w_4$, $\hat w_6$ | `op.steinhardt(l, invariant="w", normalized=True, cutoff=3.5)` |
| $\bar q$ | `op.local_average(op.tetrahedral_order(), cutoff=3.5)` |
| smooth CN | `op.smooth_coordination(r0=3.5)` |
| six-ring fraction | `op.ring_statistics(10, cutoff=3.5).fractions[:, 3]` (columns are ring sizes 3, 4, 5, …; needs hydrogens) |
| F4 | `op.F4(cutoff=3.5)` (needs hydrogens) |
| Voronoi asphericity | `op.voronoi_asphericity()` |
| CHILL, CHILL+ | `op.chill()`, `op.chill_plus()` |
| template matching | `op.template_matching()` |

Here `op = OrderParameter(u, cutoff=5.0)`, the setting the benchmark used.

## Test a choice on your own data

The quickest check: take a frame you know is ice and a frame you know is
liquid, and find the cut that separates them best. The fewer molecules end up
on the wrong side, the better the parameter.

```python
import numpy as np
import MDAnalysis as mda
from SolOrder import OrderParameter, eFingerprint

def best_threshold(ice, liquid):
    """Cut that misclassifies the fewest molecules, whichever side ice lies on."""
    best = (np.inf, None)
    for t in np.unique(np.r_[ice, liquid]):
        above = (np.mean(ice < t) + np.mean(liquid >= t)) / 2   # ice above the cut
        below = (np.mean(ice >= t) + np.mean(liquid < t)) / 2   # ice below the cut
        best = min(best, (min(above, below), t))
    return best[1], best[0]

u_ice, u_liquid = mda.Universe("ice_Ih.gro"), mda.Universe("liquid.gro")
ice, liquid = OrderParameter(u_ice, cutoff=5.0), OrderParameter(u_liquid, cutoff=5.0)
candidates = {
    "coordination number": lambda op, u: op.coordination_number(3.5),
    "q": lambda op, u: op.tetrahedral_order(),
    "averaged q": lambda op, u: op.local_average(op.tetrahedral_order(), cutoff=3.5),
    "averaged q6": lambda op, u: op.steinhardt(6, averaged=True, cutoff=3.5),
    "local entropy": lambda op, u: eFingerprint(u, cutoff=5.0).localEntropy(),
}
for name, compute in candidates.items():
    t, error = best_threshold(compute(ice, u_ice), compute(liquid, u_liquid))
    print(f"{name:20s} cut at {t:7.3f}   misclassified {100*error:4.1f}%")
```

```text
coordination number  cut at   5.000   misclassified 37.1%
q                    cut at   0.915   misclassified 15.9%
averaged q           cut at   0.919   misclassified  2.2%
averaged q6          cut at   0.389   misclassified  0.9%
local entropy        cut at  -3.513   misclassified  0.0%
```

Averaging over neighbors matters: $\bar q$ misclassifies about seven times
fewer molecules than $q$.

!!! warning "A threshold fitted on one frame is optimistic"
    This example fits and scores the cut on the same molecules. For a number
    you can report, fit the cut on some frames (or one simulation) and measure
    the error on others. Thermal motion, temperature, pressure and the water
    model all move the distributions, so a threshold rarely transfers unchanged.

!!! tip "Interfaces are harder than bulk"
    A parameter that separates bulk ice from bulk liquid perfectly can still
    mislabel the molecules at an ice surface, which are neither. If you study
    nuclei or interfaces, test on an interface too (see
    [Classifying ice](ice-classification.md)).

Next: [working with trajectories](trajectories.md).
