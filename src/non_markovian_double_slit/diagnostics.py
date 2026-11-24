"""Numerical diagnostics shared by the demonstrations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.signal import find_peaks


def l1_distance(first: ArrayLike, second: ArrayLike) -> float:
    """L1 distance between equal-shape arrays."""
    a = np.asarray(first, dtype=float)
    b = np.asarray(second, dtype=float)
    if a.shape != b.shape:
        raise ValueError("arrays must have equal shape")
    return float(np.sum(np.abs(a - b)))


def total_variation_distance(first: ArrayLike, second: ArrayLike) -> float:
    return 0.5 * l1_distance(first, second)


def periodic_bin_masses(
    centers: ArrayLike, masses: ArrayLike, detector_centers: ArrayLike
) -> np.ndarray:
    """Integrate uniform density within each input cell over new periodic bins.

    Both grids have uniform spacing and span one period. A detector bin can
    straddle the seam even when the two grids have different center offsets.
    """
    x = np.asarray(centers, dtype=float)
    p = np.asarray(masses, dtype=float)
    y = np.asarray(detector_centers, dtype=float)
    if x.ndim != 1 or p.shape != x.shape or y.ndim != 1 or min(x.size, y.size) < 2:
        raise ValueError("expected two grids and one mass per source bin")
    dx, dy = x[1] - x[0], y[1] - y[0]
    length = x.size * dx
    if not (
        dx > 0
        and dy > 0
        and np.allclose(np.diff(x), dx)
        and np.allclose(np.diff(y), dy)
        and np.isclose(length, y.size * dy)
    ):
        raise ValueError("grids must be uniform and have the same period")
    # Extend the cumulative mass by one total per complete turn.
    edges = np.r_[y - dy / 2, y[-1] + dy / 2] - (x[0] - dx / 2)
    turns, offset = np.divmod(edges, length)
    cumulative = np.r_[0.0, np.cumsum(p)]
    values = np.interp(offset, np.arange(x.size + 1) * dx, cumulative)
    return np.diff(values + turns * cumulative[-1])


def central_fringe_spacing(
    positions: ArrayLike,
    intensity: ArrayLike,
    *,
    expected_spacing: float,
    half_window: float | None = None,
) -> float:
    """Estimate median central peak separation with an expected-scale prior."""
    x = np.asarray(positions, dtype=float)
    y = np.asarray(intensity, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or x.size < 5:
        raise ValueError("positions and intensity must be equal one-dimensional arrays")
    if expected_spacing <= 0:
        raise ValueError("expected_spacing must be positive")
    if not np.all(np.diff(x) > 0):
        raise ValueError("positions must be strictly increasing")
    if half_window is None:
        # Restrict the estimate to peaks near the center of the envelope.
        half_window = 2 * expected_spacing
    mask = np.abs(x) <= half_window
    local_x = x[mask]
    local_y = y[mask]
    grid_spacing = float(np.median(np.diff(local_x)))
    minimum_peak_distance = max(1, round(0.55 * expected_spacing / grid_spacing))
    peaks, _ = find_peaks(
        local_y,
        distance=minimum_peak_distance,
        prominence=max(0.01 * float(np.max(local_y)), np.finfo(float).eps),
    )
    if peaks.size < 3:
        raise RuntimeError("fewer than three central fringes were resolved")
    return float(np.median(np.diff(local_x[peaks])))


def normalized_correlation(first: ArrayLike, second: ArrayLike) -> float:
    """Cosine similarity after removing each array's mean."""
    a = np.asarray(first, dtype=float)
    b = np.asarray(second, dtype=float)
    if a.shape != b.shape:
        raise ValueError("arrays must have equal shape")
    a = a - a.mean()
    b = b - b.mean()
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    if denominator == 0:
        raise ValueError("correlation is undefined for a constant array")
    return float(np.dot(a.ravel(), b.ravel()) / denominator)
