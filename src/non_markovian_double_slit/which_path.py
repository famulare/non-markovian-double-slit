"""Which-path records, physical erasure, and screen marginals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


def screen_pattern(
    upper_amplitude: ArrayLike,
    lower_amplitude: ArrayLike,
    *,
    upper_weight: float = 0.5,
    relative_phase: float = 0.0,
    record_overlap: complex = 1.0,
) -> NDArray[np.float64]:
    """Electron marginal after tracing out a path record.

    ``record_overlap`` is ``<r_upper|r_lower>``. Its modulus cannot exceed one.
    A value of one gives full coherence and zero gives a perfectly retained
    which-path record.
    """
    upper = np.asarray(upper_amplitude, dtype=complex)
    lower = np.asarray(lower_amplitude, dtype=complex)
    if upper.shape != lower.shape:
        raise ValueError("upper and lower amplitudes must have the same shape")
    if not 0 <= upper_weight <= 1:
        raise ValueError("upper_weight must lie in [0, 1]")
    if abs(record_overlap) > 1 + 1e-13:
        raise ValueError("record overlap cannot exceed one in modulus")
    lower_weight = 1 - upper_weight
    mixture = upper_weight * np.abs(upper) ** 2 + lower_weight * np.abs(lower) ** 2
    cross = (
        2
        * np.sqrt(upper_weight * lower_weight)
        * np.real(np.exp(1j * relative_phase) * np.conj(upper) * lower * record_overlap)
    )
    pattern = mixture + cross
    pattern[np.abs(pattern) < 1e-15] = 0.0
    if np.min(pattern) < -1e-12:
        raise RuntimeError("computed marginal is negative")
    return pattern


def retained_record_pattern(
    upper_amplitude: ArrayLike,
    lower_amplitude: ArrayLike,
    *,
    upper_weight: float = 0.5,
) -> NDArray[np.float64]:
    return screen_pattern(
        upper_amplitude,
        lower_amplitude,
        upper_weight=upper_weight,
        record_overlap=0.0,
    )


def local_reset_pattern(
    upper_amplitude: ArrayLike,
    lower_amplitude: ArrayLike,
    *,
    upper_weight: float = 0.5,
) -> NDArray[np.float64]:
    """Pattern after a physical local reset moves the record to a bath.

    A local completely positive trace-preserving reset does not destroy the
    record in the total environment. The branch overlap therefore remains zero,
    and the electron marginal is unchanged from a retained record.
    """
    return retained_record_pattern(
        upper_amplitude, lower_amplitude, upper_weight=upper_weight
    )


def coherent_uncomputation_pattern(
    upper_amplitude: ArrayLike,
    lower_amplitude: ArrayLike,
    *,
    upper_weight: float = 0.5,
    relative_phase: float = 0.0,
) -> NDArray[np.float64]:
    """Pattern after reversing the path-record interaction coherently."""
    return screen_pattern(
        upper_amplitude,
        lower_amplitude,
        upper_weight=upper_weight,
        relative_phase=relative_phase,
        record_overlap=1.0,
    )


@dataclass(frozen=True)
class EraserPatterns:
    """Joint screen probabilities for the ``+`` and ``-`` record outcomes."""

    plus: NDArray[np.float64]
    minus: NDArray[np.float64]

    @property
    def aggregate(self) -> NDArray[np.float64]:
        return self.plus + self.minus


def quantum_eraser_patterns(
    upper_amplitude: ArrayLike,
    lower_amplitude: ArrayLike,
    *,
    upper_weight: float = 0.5,
    relative_phase: float = 0.0,
) -> EraserPatterns:
    """Unnormalized joint patterns conditioned on the ``+/-`` record basis."""
    upper = np.asarray(upper_amplitude, dtype=complex)
    lower = np.asarray(lower_amplitude, dtype=complex)
    if upper.shape != lower.shape:
        raise ValueError("upper and lower amplitudes must have the same shape")
    branch_upper = np.sqrt(upper_weight) * upper
    branch_lower = np.exp(1j * relative_phase) * np.sqrt(1 - upper_weight) * lower
    plus = 0.5 * np.abs(branch_upper + branch_lower) ** 2
    minus = 0.5 * np.abs(branch_upper - branch_lower) ** 2
    return EraserPatterns(plus=plus, minus=minus)


def path_distinguishability(record_overlap: complex) -> float:
    """Pure-record distinguishability ``sqrt(1 - |overlap|^2)``."""
    if abs(record_overlap) > 1 + 1e-13:
        raise ValueError("record overlap cannot exceed one in modulus")
    return float(np.sqrt(max(0.0, 1 - abs(record_overlap) ** 2)))
