"""Hamiltonian construction of finite unistochastic double-slit laws."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .process import FiniteIndivisibleProcess

ComplexVector = NDArray[np.complex128]
ComplexMatrix = NDArray[np.complex128]


def periodic_grid(size: int, length: float) -> NDArray[np.float64]:
    """Centered cell locations for an odd periodic grid."""
    if size < 3 or size % 2 == 0:
        raise ValueError("size must be odd and at least three")
    if length <= 0:
        raise ValueError("length must be positive")
    spacing = length / size
    return (np.arange(size) - size // 2) * spacing


def normalized_box_mode(grid: ArrayLike, center: float, width: float) -> ComplexVector:
    """Normalized finite-state aperture mode supported on grid cells."""
    x = np.asarray(grid, dtype=float)
    if width <= 0:
        raise ValueError("width must be positive")
    mask = np.abs(x - center) <= width / 2
    count = int(mask.sum())
    if count == 0:
        raise ValueError("aperture contains no grid cells")
    mode = mask.astype(complex) / np.sqrt(count)
    return mode


def unitary_mapping(source: ArrayLike, target: ArrayLike) -> ComplexMatrix:
    """Return a unitary that maps one normalized vector to another.

    This complex Householder map is a convenient completion; only the target
    column is physically used for the prepared source. Building the full unitary
    nevertheless makes every column of ``Gamma = abs(U)**2`` auditable.
    """
    a = np.asarray(source, dtype=complex)
    b = np.asarray(target, dtype=complex)
    if a.ndim != 1 or b.shape != a.shape:
        raise ValueError("source and target must be vectors of the same shape")
    if not np.isclose(np.vdot(a, a), 1.0, rtol=0, atol=1e-12):
        raise ValueError("source must be normalized")
    if not np.isclose(np.vdot(b, b), 1.0, rtol=0, atol=1e-12):
        raise ValueError("target must be normalized")
    identity = np.eye(a.size, dtype=complex)
    if np.array_equal(a, b):
        return identity
    difference = a - b
    if np.linalg.norm(difference) < 1e-3:
        # Reflect a to -a first, then map -a to b with a well-scaled
        # denominator. Near coincidence the one-reflection formula loses
        # the O(||a-b||²) real part of its denominator to cancellation.
        v = a + b
        second = identity - np.outer(v, v.conj()) / np.vdot(v, a)
        return second - 2 * np.outer(second @ a, a.conj())
    return identity - np.outer(difference, difference.conj()) / np.vdot(difference, a)


def free_propagator(
    size: int,
    length: float,
    elapsed: float,
    *,
    mass: float = 1.0,
    hbar: float = 1.0,
) -> ComplexMatrix:
    """Exact spectral propagator for the finite periodic kinetic Hamiltonian."""
    if mass <= 0 or hbar <= 0:
        raise ValueError("mass and hbar must be positive")
    spacing = length / size
    wave_numbers = 2 * np.pi * np.fft.fftfreq(size, d=spacing)
    phase = np.exp(-0.5j * hbar * wave_numbers**2 * elapsed / mass)
    basis = np.eye(size, dtype=complex)
    momentum_columns = np.fft.fft(basis, axis=0, norm="ortho")
    return np.fft.ifft(phase[:, np.newaxis] * momentum_columns, axis=0, norm="ortho")


def unitarity_error(matrix: ArrayLike) -> float:
    values = np.asarray(matrix, dtype=complex)
    identity = np.eye(values.shape[0])
    return float(np.max(np.abs(values.conj().T @ values - identity)))


@dataclass(frozen=True)
class FiniteDoubleSlitModel:
    """Finite cyclic-grid double slit with one definite source configuration."""

    grid: NDArray[np.float64]
    upper_mode: ComplexVector
    lower_mode: ComplexVector
    source_index: int
    upper_weight: float
    relative_phase: float
    mass: float
    hbar: float
    preparation: ComplexMatrix
    source_time: float = -1.0
    slit_time: float = 0.0

    @classmethod
    def create(
        cls,
        *,
        size: int = 101,
        length: float = 40.0,
        slit_separation: float = 6.0,
        slit_width: float = 1.5,
        upper_weight: float = 0.5,
        relative_phase: float = 0.0,
        source_index: int | None = None,
        mass: float = 1.0,
        hbar: float = 1.0,
    ) -> FiniteDoubleSlitModel:
        if not 0 <= upper_weight <= 1:
            raise ValueError("upper_weight must lie in [0, 1]")
        if slit_separation <= slit_width:
            raise ValueError("slits must be disjoint")
        grid = periodic_grid(size, length)
        upper = normalized_box_mode(grid, slit_separation / 2, slit_width)
        lower = normalized_box_mode(grid, -slit_separation / 2, slit_width)
        if abs(np.vdot(upper, lower)) > 1e-14:
            raise ValueError("finite aperture modes overlap")
        if source_index is None:
            source_index = size // 2
        if not 0 <= source_index < size:
            raise ValueError("source_index is out of bounds")
        state = np.sqrt(upper_weight) * upper + (
            np.exp(1j * relative_phase) * np.sqrt(1 - upper_weight) * lower
        )
        source = np.zeros(size, dtype=complex)
        source[source_index] = 1.0
        preparation = unitary_mapping(source, state)
        return cls(
            grid=grid,
            upper_mode=upper,
            lower_mode=lower,
            source_index=source_index,
            upper_weight=upper_weight,
            relative_phase=relative_phase,
            mass=mass,
            hbar=hbar,
            preparation=preparation,
        )

    @property
    def size(self) -> int:
        return self.grid.size

    @property
    def length(self) -> float:
        return float(self.grid.size * (self.grid[1] - self.grid[0]))

    @property
    def source_state(self) -> ComplexVector:
        state = np.zeros(self.size, dtype=complex)
        state[self.source_index] = 1.0
        return state

    @property
    def slit_state(self) -> ComplexVector:
        return self.preparation @ self.source_state

    def propagator(self, elapsed: float) -> ComplexMatrix:
        if elapsed < 0:
            raise ValueError("elapsed time must be nonnegative")
        return free_propagator(
            self.size,
            self.length,
            elapsed,
            mass=self.mass,
            hbar=self.hbar,
        )

    def evolution(self, elapsed: float) -> ComplexMatrix:
        """Unitary from the definite source to a target after the slits."""
        return self.propagator(elapsed) @ self.preparation

    def gamma(self, elapsed: float) -> NDArray[np.float64]:
        """Full positive law, including columns not used by the chosen source."""
        return np.abs(self.evolution(elapsed)) ** 2

    def probability(self, elapsed: float) -> NDArray[np.float64]:
        return self.gamma(elapsed)[:, self.source_index]

    def slit_components(self, elapsed: float) -> tuple[ComplexVector, ComplexVector]:
        propagator = self.propagator(elapsed)
        return propagator @ self.upper_mode, propagator @ self.lower_mode

    def divided_law(
        self, target_elapsed: float, intermediate_elapsed: float
    ) -> NDArray[np.float64]:
        """Law after recording the intermediate position to single-cell resolution."""
        if not 0 <= intermediate_elapsed <= target_elapsed:
            raise ValueError("require 0 <= intermediate <= target")
        later = np.abs(self.propagator(target_elapsed - intermediate_elapsed)) ** 2
        earlier = np.abs(self.evolution(intermediate_elapsed)) ** 2
        return later @ earlier

    def interference_residual(
        self, target_elapsed: float, intermediate_elapsed: float
    ) -> NDArray[np.float64]:
        return self.gamma(target_elapsed) - self.divided_law(
            target_elapsed, intermediate_elapsed
        )

    def process(self, elapsed_times: Iterable[float]) -> FiniteIndivisibleProcess:
        """Cross the amplitude boundary and return only positive laws."""
        elapsed = tuple(sorted(set(float(t) for t in elapsed_times)))
        if any(t < 0 for t in elapsed):
            raise ValueError("elapsed times must be nonnegative")
        targets = tuple(self.slit_time + t for t in elapsed)
        times = tuple(dict.fromkeys((self.source_time, self.slit_time, *targets)))
        laws: dict[tuple[float, float], NDArray[np.float64]] = {
            (self.source_time, self.source_time): np.eye(self.size),
            (self.slit_time, self.source_time): np.abs(self.preparation) ** 2,
        }
        for target, duration in zip(targets, elapsed, strict=True):
            laws[(target, self.source_time)] = self.gamma(duration)
        initial = np.zeros(self.size)
        initial[self.source_index] = 1.0
        return FiniteIndivisibleProcess(
            configurations=tuple(f"x={x:.12g}" for x in self.grid),
            times=times,
            conditioning_times=frozenset({self.source_time}),
            laws=laws,
            initial_time=self.source_time,
            initial_distribution=initial,
        )
