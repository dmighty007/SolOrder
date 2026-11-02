"""Numba kernels. All inputs are numerical snapshots; no Universe is accessed."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def coordination_chunk(start, stop, offsets, distances, radius, smooth, n, m):
    out = np.zeros(stop - start)
    for i in range(start, stop):
        value = 0.
        for p in range(offsets[i], offsets[i + 1]):
            x = distances[p] / radius
            if not smooth:
                value += x <= 1.
            elif abs(x - 1.) < 1e-12:
                value += n / m
            elif x == 0.:
                value += 1.
            elif x < 1.:
                value += np.expm1(n * np.log(x)) / np.expm1(m * np.log(x))
            else:
                # Equivalent inverse-power form avoids overflowing at large x.
                value += x**(n-m) * (-np.expm1(-n*np.log(x))) / (-np.expm1(-m*np.log(x)))
        out[i - start] = value
    return out


@njit(cache=True, nogil=True)
def kth_chunk(start, stop, offsets, distances, k):
    out = np.full(stop - start, np.nan)
    for i in range(start, stop):
        if offsets[i + 1] - offsets[i] >= k:
            out[i - start] = distances[offsets[i] + k - 1]
    return out


@njit(cache=True, nogil=True)
def tetrahedral_chunk(start, stop, offsets, distances, vectors, orientational):
    out = np.full(stop - start, np.nan)
    for i in range(start, stop):
        p = offsets[i]
        if offsets[i + 1] - p < 4 or distances[p] <= 0:
            continue
        if orientational:
            total = 0.
            for j in range(3):
                for k in range(j + 1, 4):
                    cosine = np.dot(vectors[p+j], vectors[p+k]) / (distances[p+j] * distances[p+k])
                    total += (min(1., max(-1., cosine)) + 1./3.)**2
            out[i - start] = 1. - 3./8. * total
        else:
            radii = distances[p:p+4]
            mean = np.mean(radii)
            out[i - start] = 1. - np.sum((radii - mean)**2) / (12. * mean**2)
    return out


@njit(cache=True, nogil=True)
def average_chunk(start, stop, offsets, indices, values):
    out = np.empty((stop - start, values.shape[1]), dtype=values.dtype)
    for i in range(start, stop):
        row = values[i].copy()
        for p in range(offsets[i], offsets[i + 1]):
            row += values[indices[p]]
        out[i - start] = row / (offsets[i + 1] - offsets[i] + 1)
    return out


@njit(cache=True, nogil=True)
def lsi_chunk(start, stop, offsets, distances, shell):
    out = np.full(stop - start, np.nan)
    for i in range(start, stop):
        radii = distances[offsets[i]:offsets[i+1]]
        n = 0
        while n < len(radii) and radii[n] < shell:
            n += 1
        # The first distance outside the shell is needed for the final gap.
        if n > 0 and n < len(radii):
            gaps = radii[1:n+1] - radii[:n]
            out[i-start] = np.mean((gaps - np.mean(gaps))**2)
    return out


@njit(cache=True, nogil=True)
def harmonic_mean(offsets, harmonics):
    out = np.full((len(offsets)-1, harmonics.shape[1]), np.nan + 0j)
    for i in range(len(out)):
        if offsets[i+1] > offsets[i]:
            out[i] = np.sum(harmonics[offsets[i]:offsets[i+1]], axis=0) / (offsets[i+1]-offsets[i])
    return out


@njit(cache=True, nogil=True)
def third_order_chunk(start, stop, qlm, terms, coefficients, normalized):
    out = np.full(stop-start, np.nan)
    for i in range(start, stop):
        norm = np.sum(np.abs(qlm[i])**2)
        if norm == 0 or not np.isfinite(norm):
            continue
        value = 0j
        for p in range(len(coefficients)):
            a, b, c = terms[p]
            value += coefficients[p] * qlm[i, a] * qlm[i, b] * qlm[i, c]
        out[i-start] = value.real / norm**1.5 if normalized else value.real
    return out


@njit(cache=True, nogil=True)
def connections_chunk(start, stop, offsets, indices, qlm, threshold):
    out = np.zeros(stop-start, dtype=np.int64)
    for i in range(start, stop):
        ni = np.sqrt(np.sum(np.abs(qlm[i])**2))
        if ni == 0 or not np.isfinite(ni):
            continue
        for p in range(offsets[i], offsets[i+1]):
            j = indices[p]
            nj = np.sqrt(np.sum(np.abs(qlm[j])**2))
            if nj > 0 and np.isfinite(nj):
                correlation = np.sum(qlm[i] * np.conj(qlm[j])).real / (ni * nj)
                out[i-start] += correlation > threshold
    return out


@njit(cache=True, nogil=True)
def hbond_chunk(start, stop, source, indices, vectors, distances, oh, cutoff, cosine):
    out = np.zeros((stop-start, 2), dtype=np.bool_)
    for p in range(start, stop):
        if distances[p] == 0 or distances[p] > cutoff:
            continue
        i = source[p]
        for h in range(2):
            alignment = np.dot(oh[i, h], vectors[p]) / (np.linalg.norm(oh[i, h]) * distances[p])
            out[p-start, h] = alignment > cosine
    return out


@njit(cache=True, nogil=True)
def dihedral_value(a, b, c, d):
    n1 = np.cross(b-a, c-b)
    n2 = np.cross(c-b, d-c)
    denom = np.linalg.norm(n1) * np.linalg.norm(n2)
    if denom == 0:
        return np.nan
    return -np.arccos(min(1., max(-1., np.dot(n1, n2) / denom)))


@njit(cache=True, nogil=True)
def water_dihedral(oo, ha, hb):
    # Choose the farther hydrogen independently on each molecule. A global
    # minimum/maximum among all four distances can select the wrong hydrogen.
    ia = 0 if np.linalg.norm(ha[0]-oo) >= np.linalg.norm(ha[1]-oo) else 1
    ib = 0 if np.linalg.norm(oo+hb[0]) >= np.linalg.norm(oo+hb[1]) else 1
    return dihedral_value(ha[ia], np.zeros(3), oo, oo+hb[ib])


@njit(cache=True, nogil=True)
def dihedrals_chunk(start, stop, source, indices, vectors, oh):
    out = np.empty(stop-start)
    for p in range(start, stop):
        out[p-start] = water_dihedral(vectors[p], oh[source[p]], oh[indices[p]])
    return out


@njit(cache=True, nogil=True)
def row_mean_chunk(start, stop, offsets, values):
    out = np.full(stop-start, np.nan)
    for i in range(start, stop):
        if offsets[i+1] > offsets[i]:
            out[i-start] = np.mean(values[offsets[i]:offsets[i+1]])
    return out


@njit(cache=True, nogil=True)
def minimum_angle_chunk(start, stop, offsets, indices, vectors, distances, oh):
    out = np.full(stop-start, np.nan)
    for i in range(start, stop):
        p = offsets[i]
        if p == offsets[i+1] or distances[p] == 0:
            continue
        j = indices[p]
        best = -1.
        for h in range(2):
            best = max(best, np.dot(oh[i, h], vectors[p]) / (np.linalg.norm(oh[i, h])*distances[p]))
            best = max(best, -np.dot(oh[j, h], vectors[p]) / (np.linalg.norm(oh[j, h])*distances[p]))
        out[i-start] = np.arccos(min(1., max(-1., best)))
    return out


@njit(cache=True, nogil=True)
def entropy_chunk(start, stop, centers, offsets, distances, r, sigma, global_rho, local, cutoff):
    out = np.full(stop-start, np.nan)
    for row in range(start, stop):
        i = centers[row]
        radii = distances[offsets[i]:offsets[i+1]]
        rho = len(radii) / (4./3.*np.pi*cutoff**3) if local else global_rho
        if rho == 0:
            continue
        integrand = np.zeros(len(r))
        for k in range(1, len(r)):
            g = 0.
            for rij in radii:
                g += np.exp(-(r[k]-rij)**2 / (2.*sigma**2))
            g /= 4.*np.pi*rho*r[k]**2 * np.sqrt(2.*np.pi*sigma**2)
            integrand[k] = (g*np.log(g)-g+1.)*r[k]**2 if g > 1e-10 else r[k]**2
        integral = 0.
        for k in range(1, len(r)):
            integral += .5*(integrand[k-1]+integrand[k])*(r[k]-r[k-1])
        out[row-start] = -2.*np.pi*rho*integral
    return out
