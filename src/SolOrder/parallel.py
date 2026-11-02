"""Ordered joblib execution of array kernels, without nested Numba thread pools."""

from numbers import Integral

import numpy as np
from joblib import Parallel, delayed, effective_n_jobs


def validate_parallel(n_jobs, backend):
    if isinstance(n_jobs, bool) or not isinstance(n_jobs, Integral) or n_jobs == 0:
        raise ValueError("n_jobs must be a nonzero integer")
    if backend not in ("threading", "loky"):
        raise ValueError("backend must be 'threading' or 'loky'")


def map_chunks(kernel, size, *args, n_jobs=1, backend="threading"):
    """Concatenate kernel(start, stop, *args) results in input order.

    Kernels must return an array with one row per input item. The first chunk
    runs on the caller, compiling Numba kernels before threads start. Numerical
    kernels release the GIL and do not start their own parallel thread pools.
    ``loky`` is available for process isolation; threads share arrays directly.
    """
    validate_parallel(n_jobs, backend)
    if size == 0:
        return kernel(0, 0, *args)
    workers = min(size, effective_n_jobs(n_jobs))
    if workers == 1:
        return kernel(0, size, *args)
    edges = np.linspace(0, size, min(size, 4 * workers) + 1, dtype=int)
    first = kernel(int(edges[0]), int(edges[1]), *args)
    rest = Parallel(n_jobs=workers, backend=backend)(
        delayed(kernel)(int(start), int(stop), *args)
        for start, stop in zip(edges[1:-1], edges[2:])
    )
    return np.concatenate([first, *rest], axis=0)
