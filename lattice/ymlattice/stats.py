"""Error analysis helpers: binning, jackknife, integrated autocorrelation time."""

from __future__ import annotations

import numpy as np


def binned_mean_err(x: np.ndarray, nbins: int = 20):
    """Mean and standard error from ``nbins`` consecutive bins (robust against autocorrelation)."""
    x = np.asarray(x, dtype=float)
    n = len(x) // nbins * nbins
    if n < nbins:
        return float(np.mean(x)), float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else np.nan
    b = x[:n].reshape(nbins, -1).mean(axis=1)
    return float(b.mean()), float(b.std(ddof=1) / np.sqrt(nbins))


def jackknife(samples: np.ndarray, func, nblocks: int = 20):
    """Jackknife estimate of ``func(samples)`` with blocks of consecutive samples.

    ``samples`` has shape ``(nmeas, ...)``; ``func`` maps an array of samples to a (vector of)
    secondary quantities.  Returns ``(estimate, error)``."""
    samples = np.asarray(samples)
    n = len(samples) // nblocks * nblocks
    samples = samples[:n]
    blocks = samples.reshape(nblocks, -1, *samples.shape[1:])
    full = func(samples)
    ests = []
    for i in range(nblocks):
        keep = np.concatenate([blocks[:i], blocks[i + 1:]], axis=0).reshape(-1, *samples.shape[1:])
        ests.append(func(keep))
    ests = np.array(ests)
    err = np.sqrt((nblocks - 1) / nblocks * np.sum((ests - ests.mean(axis=0)) ** 2, axis=0))
    return full, err


def autocorrelation_time(x: np.ndarray, wmax: int | None = None) -> float:
    """Integrated autocorrelation time ``tau_int = 1/2 + sum_{t>=1} rho(t)`` with automatic windowing."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    var = np.dot(x, x) / n
    if var == 0:
        return 0.5
    wmax = wmax or n // 4
    tau = 0.5
    for t in range(1, wmax):
        rho = np.dot(x[:-t], x[t:]) / ((n - t) * var)
        if rho <= 0:
            break
        tau += rho
        if t >= 6 * tau:  # automatic window (Sokal)
            break
    return float(tau)
