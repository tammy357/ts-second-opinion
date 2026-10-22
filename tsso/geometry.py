"""Geometry helpers: distance-matrix error (D-MAE), model disagreement, consensus pick."""
import numpy as np


def pair_distances(xyz):
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    return d[np.triu_indices(len(xyz), 1)]


def dmae(a, b):
    """Mean absolute error between the interatomic distances of two geometries (GoFlow's D-MAE)."""
    return float(np.mean(np.abs(pair_distances(a) - pair_distances(b))))


def disagreement_matrix(preds):
    D = np.array([pair_distances(p) for p in preds])
    return np.abs(D[:, None, :] - D[None, :, :]).mean(-1)


def consensus_index(preds, allowed=None):
    """Index of the prediction closest to the element-wise median of all predictions' distances.
    If `allowed` (bool mask) is given, choose only among allowed predictions (falls back to all)."""
    D = np.array([pair_distances(p) for p in preds])
    med = np.median(D, axis=0)
    score = np.abs(D - med).mean(1)
    if allowed is not None and np.any(allowed):
        score = np.where(allowed, score, np.inf)
    return int(np.argmin(score))


def disputed_pair(preds):
    """Atom pair whose distance varies most across predictions, with its spread (Angstrom)."""
    n = preds.shape[1]
    d = np.linalg.norm(preds[:, :, None, :] - preds[:, None, :, :], axis=-1)
    sd = d.std(0)
    np.fill_diagonal(sd, 0)
    i, j = np.unravel_index(np.argmax(sd), sd.shape)
    return (int(min(i, j)), int(max(i, j))), float(sd[i, j]), d[:, i, j]


def signed_volume(xyz, centre, nbrs):
    v = xyz[nbrs[:3]] - xyz[centre]
    return float(np.dot(v[0], np.cross(v[1], v[2])))
