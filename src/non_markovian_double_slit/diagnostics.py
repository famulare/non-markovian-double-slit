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
        # Use the innermost fringes: farther peaks acquire the expected
        # nonparaxial angular-spacing correction.
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
