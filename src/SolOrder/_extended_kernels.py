"""Numerical kernels for angular, symmetry and ice-classification measures."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def centrosymmetry_chunk(start, stop, offsets, vectors, count, mode):
    out = np.full(stop-start, np.nan)
    for i in range(start, stop):
        first = offsets[i]
        if offsets[i+1]-first < count:
            continue
        v = vectors[first:first+count]
        if mode == 1:  # Manuscript's distance-rank pairing.
            value = 0.
            for a in range(count//2):
                value += np.sum((v[a]+v[a+count//2])**2)
        else:  # Conventional Kelchner/LAMMPS definition.
            weights = np.empty(count*(count-1)//2)
            p = 0
            for a in range(count):
                for b in range(a+1, count):
                    weights[p] = np.sum((v[a]+v[b])**2)
                    p += 1
            value = np.sum(np.sort(weights)[:count//2])
        out[i-start] = value
    return out


@njit(cache=True, nogil=True)
def wedge_chunk(start, stop, left, right, vectors, distances):
    out = np.empty((stop-start, 3))
    for p in range(start, stop):
        a, b = left[p], right[p]
        ra, rb = distances[a], distances[b]
        out[p-start, 0] = ra
        out[p-start, 1] = rb
        out[p-start, 2] = np.arccos(min(1., max(-1., np.dot(vectors[a], vectors[b])/(ra*rb))))
    return out


@njit(cache=True, nogil=True)
def pair_angle_chunk(start, stop, source, indices, vectors, distances, oh):
    out = np.full(stop-start, np.nan)
    for p in range(start, stop):
        if distances[p] <= 0:
            continue
        i, j = source[p], indices[p]
        best = -1.
        for h in range(2):
            best = max(best, np.dot(oh[i, h], vectors[p])/(np.linalg.norm(oh[i, h])*distances[p]))
            best = max(best, -np.dot(oh[j, h], vectors[p])/(np.linalg.norm(oh[j, h])*distances[p]))
        out[p-start] = np.arccos(min(1., max(-1., best)))
    return out


@njit(cache=True, nogil=True)
def classify_chunk(start, stop, offsets, indices, staggered, eclipsed, plus):
    # 0 other/liquid, 1 Ih, 2 Ic, 3 interfacial ice, 4 hydrate, 5 interface hydrate.
    out = np.zeros(stop-start, dtype=np.int64)
    for i in range(start, stop):
        if offsets[i+1]-offsets[i] != 4:
            continue
        s, e = staggered[i], eclipsed[i]
        if s == 4:
            out[i-start] = 2
        elif s == 3 and e == 1:
            out[i-start] = 1
        elif plus:
            if s == 2 or s == 3:
                out[i-start] = 3
            elif e == 4:
                out[i-start] = 4
            elif e == 3:
                out[i-start] = 5
        else:
            for p in range(offsets[i], offsets[i+1]):
                sj = staggered[indices[p]]
                if (s == 2 and sj > 2) or (s == 3 and e == 0 and sj == 2):
                    out[i-start] = 3
                    break
    return out


@njit(cache=True, nogil=True)
def minimum_image(d, cell, inverse, half_height):
    """Shortest periodic image of d (cell rows are box vectors).

    A wrapped vector shorter than half the smallest cell height is already the
    minimum image; otherwise the 27 neighboring images are searched.
    """
    frac = d @ inverse
    frac -= np.round(frac)
    base = frac @ cell
    best_sq = np.sum(base**2)
    if best_sq < half_height**2:
        return base
    best = base.copy()
    for a in range(-1, 2):
        for b in range(-1, 2):
            for c in range(-1, 2):
                trial = base + a*cell[0] + b*cell[1] + c*cell[2]
                sq = np.sum(trial**2)
                if sq < best_sq:
                    best, best_sq = trial, sq
    return best


@njit(cache=True, nogil=True)
def neighborhood_chunk(start, stop, offsets, vectors, distances, cell, inverse, half_height, m):
    out = np.full((stop-start, 3), np.nan)
    for i in range(start, stop):
        lo, hi = offsets[i], offsets[i+1]
        n = hi-lo
        if n < 2:
            continue
        v = vectors[lo:hi]
        delta = np.empty((n, n, 3))
        score = np.empty((n, n))
        for j in range(n):
            delta[j, j] = 0.
            score[j, j] = np.inf
            for k in range(j+1, n):
                d = minimum_image(v[k]-v[j], cell, inverse, half_height)
                delta[j, k] = d
                delta[k, j] = -d
                r = np.sqrt(np.sum(d**2))
                score[j, k] = distances[lo+k]+r
                score[k, j] = distances[lo+j]+r
        ranks = np.empty((n, n), dtype=np.int64)
        for j in range(n):
            ranks[j] = np.argsort(score[j], kind="mergesort")
        out[i-start] = neighborhood_reduce(v, delta, ranks, m)
    return out


@njit(cache=True, nogil=True)
def neighborhood_reduce(vectors, delta, ranks, m):
    n = len(vectors)
    out = np.zeros(3)
    total = np.zeros(3)
    for j in range(n):
        sum_a = np.zeros(3)
        sum_p = np.zeros(3)
        for k in ranks[j, :min(m, n-1)]:
            sum_a += vectors[k]+delta[j, k]
            sum_p += vectors[j]-delta[j, k]
        out[0] += np.sum(sum_a**2)
        out[1] += np.sum(sum_p**2)
        total += sum_p
    out[2] = np.sum(total**2)
    return out/n


@njit(cache=True, nogil=True)
def angle_histogram_chunk(start, stop, offsets, angles, bins, frequencies):
    out = np.full((stop-start, len(frequencies)), np.nan)
    for i in range(start, stop):
        count = offsets[i+1]-offsets[i]
        if count == 0:
            continue
        hist = np.zeros(bins)
        for p in range(offsets[i], offsets[i+1]):
            index = min(bins-1, int(angles[p]/np.pi*bins))
            hist[index] += 1.
        for f in range(len(frequencies)):
            real, imag = 0., 0.
            for b in range(bins):
                phase = -2.*np.pi*frequencies[f]*b/bins
                real += hist[b]*np.cos(phase)
                imag += hist[b]*np.sin(phase)
            out[i-start, f] = np.sqrt(real**2+imag**2)/count
    return out


@njit(cache=True, nogil=True)
def edge_correlation_chunk(start, stop, source, indices, qlm):
    out = np.full(stop-start, np.nan)
    for p in range(start, stop):
        i, j = source[p], indices[p]
        ni = np.sqrt(np.sum(np.abs(qlm[i])**2))
        nj = np.sqrt(np.sum(np.abs(qlm[j])**2))
        if ni > 0 and nj > 0 and np.isfinite(ni) and np.isfinite(nj):
            out[p-start] = np.sum(qlm[i]*np.conj(qlm[j])).real/(ni*nj)
    return out


@njit(cache=True, nogil=True)
def smallest_ring_chunk(start, stop, offsets, indices, vectors, max_size):
    """Smallest ring through each angle j-i-k, by breadth-first search j -> k avoiding i.

    Column s counts angles whose smallest ring has s members; column 0 counts
    angles with no ring up to max_size. Paths carry summed bond vectors, so a
    node reached through another periodic image is a distinct node and a path
    that wraps around the box never closes a ring.
    """
    n = len(offsets)-1
    slots = 4  # periodic images remembered per node
    out = np.zeros((stop-start, max_size+1), dtype=np.int64)
    seen = np.zeros(n, dtype=np.int64)
    images = np.empty((n, slots, 3))
    touched = np.empty(n, dtype=np.int64)
    queue_node = np.empty(n*slots, dtype=np.int64)
    queue_depth = np.empty(n*slots, dtype=np.int64)
    queue_vector = np.empty((n*slots, 3))
    for i in range(start, stop):
        for p in range(offsets[i], offsets[i+1]):
            for q in range(p+1, offsets[i+1]):
                j, k = indices[p], indices[q]
                seen[j] = 1
                images[j, 0] = vectors[p]
                touched[0] = j
                n_touched = 1
                queue_node[0], queue_depth[0], queue_vector[0] = j, 0, vectors[p]
                head, tail, found = 0, 1, -1
                while head < tail and found < 0:
                    u, depth, du = queue_node[head], queue_depth[head], queue_vector[head]
                    head += 1
                    if depth >= max_size-2:
                        continue
                    for r in range(offsets[u], offsets[u+1]):
                        v = indices[r]
                        if v == i:
                            continue
                        dv = du + vectors[r]
                        if v == k and np.sum((dv-vectors[q])**2) < 1.:
                            found = depth+1
                            break
                        new = seen[v] < slots
                        for s in range(seen[v]):
                            if np.sum((images[v, s]-dv)**2) < 1.:
                                new = False
                                break
                        if not new:
                            continue
                        if seen[v] == 0:
                            touched[n_touched] = v
                            n_touched += 1
                        images[v, seen[v]] = dv
                        seen[v] += 1
                        queue_node[tail], queue_depth[tail], queue_vector[tail] = v, depth+1, dv
                        tail += 1
                for t in range(n_touched):
                    seen[touched[t]] = 0
                out[i-start, found+2 if found >= 0 else 0] += 1
    return out
