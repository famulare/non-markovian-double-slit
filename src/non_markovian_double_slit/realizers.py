"""Positive path-space realizers for prescribed one-time process marginals."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .process import FiniteIndivisibleProcess


def monotone_transport_coupling(
    source: ArrayLike, target: ArrayLike
) -> NDArray[np.float64]:
    """One-dimensional quadratic-cost optimal coupling.

    Rows index the target state and columns the source state. For ordered
    one-dimensional support, the monotone coupling minimizes every convex
    function of displacement, including squared displacement.
    """
    source_values = np.asarray(source, dtype=float)
    target_values = np.asarray(target, dtype=float)
    if source_values.shape != target_values.shape or source_values.ndim != 1:
        raise ValueError("source and target must be vectors of equal length")
    if (
        np.min(source_values) < -1e-13
        or np.min(target_values) < -1e-13
        or not np.isclose(source_values.sum(), 1.0, atol=1e-11)
        or not np.isclose(target_values.sum(), 1.0, atol=1e-11)
    ):
        raise ValueError("source and target must be probability vectors")

    coupling = np.zeros((source_values.size, source_values.size))
    remaining_source = source_values.copy()
    remaining_target = target_values.copy()
    i = j = 0
    tolerance = 1e-15
    while i < source_values.size and j < source_values.size:
        amount = min(remaining_source[j], remaining_target[i])
        coupling[i, j] += amount
        remaining_source[j] -= amount
        remaining_target[i] -= amount
        if remaining_source[j] <= tolerance:
            j += 1
        if remaining_target[i] <= tolerance:
            i += 1
    coupling[np.abs(coupling) < tolerance] = 0.0
    if not np.allclose(coupling.sum(axis=0), source_values, atol=1e-12):
        raise RuntimeError("transport coupling has incorrect source marginal")
    if not np.allclose(coupling.sum(axis=1), target_values, atol=1e-12):
        raise RuntimeError("transport coupling has incorrect target marginal")
    return coupling


class PathRealizer(Protocol):
    process: FiniteIndivisibleProcess
    target_times: tuple[float, ...]

    def sample(self, count: int, *, rng: np.random.Generator) -> NDArray[np.int64]: ...

    def expected_squared_step(self, positions: ArrayLike) -> float: ...


def _draw(probabilities: NDArray[np.float64], rng: np.random.Generator) -> int:
    return int(rng.choice(probabilities.size, p=probabilities))


@dataclass(frozen=True)
class IndependentPathRealizer:
    """Target configurations independent conditional on the initial state."""

    process: FiniteIndivisibleProcess
    target_times: tuple[float, ...]

    def __post_init__(self) -> None:
        for target in self.target_times:
            self.process.law(target)

    def sample(self, count: int, *, rng: np.random.Generator) -> NDArray[np.int64]:
        if count < 0:
            raise ValueError("count must be nonnegative")
        paths = np.empty((count, 1 + len(self.target_times)), dtype=np.int64)
        for row in range(count):
            initial = _draw(self.process.initial_distribution, rng)
            paths[row, 0] = initial
            for column, target in enumerate(self.target_times, start=1):
                paths[row, column] = _draw(self.process.law(target)[:, initial], rng)
        return paths

    def expected_squared_step(self, positions: ArrayLike) -> float:
        x = np.asarray(positions, dtype=float)
        if x.shape != (self.process.size,):
            raise ValueError("positions has the wrong shape")
        total = 0.0
        for initial, weight in enumerate(self.process.initial_distribution):
            if weight == 0:
                continue
            previous = np.zeros(self.process.size)
            previous[initial] = 1.0
            for target in self.target_times:
                current = self.process.law(target)[:, initial]
                joint = np.outer(current, previous)
                total += weight * float(np.sum(joint * (x[:, None] - x[None, :]) ** 2))
                previous = current
        return total


@dataclass(frozen=True)
class TransportPathRealizer:
    """Conditionally comonotone, spatially smooth path coupling.

    Conditional on the initial configuration, one shared uniform rank is mapped
    through every target marginal's quantile function. This is projectively
    consistent when target times are added or removed, and every pairwise
    coupling is the one-dimensional monotone optimal transport coupling.
    """

    process: FiniteIndivisibleProcess
    target_times: tuple[float, ...]

    def __post_init__(self) -> None:
        for target in self.target_times:
            self.process.law(target)

    def _couplings(self, initial: int) -> tuple[NDArray[np.float64], ...]:
        previous = np.zeros(self.process.size)
        previous[initial] = 1.0
        result = []
        for target in self.target_times:
            current = self.process.law(target)[:, initial]
            result.append(monotone_transport_coupling(previous, current))
            previous = current
        return tuple(result)

    def sample(self, count: int, *, rng: np.random.Generator) -> NDArray[np.int64]:
        if count < 0:
            raise ValueError("count must be nonnegative")
        paths = np.empty((count, 1 + len(self.target_times)), dtype=np.int64)
        for row in range(count):
            initial = _draw(self.process.initial_distribution, rng)
            paths[row, 0] = initial
            rank = rng.random()
            for column, target in enumerate(self.target_times, start=1):
                cumulative = np.cumsum(self.process.law(target)[:, initial])
                paths[row, column] = min(
                    int(np.searchsorted(cumulative, rank, side="right")),
                    self.process.size - 1,
                )
        return paths

    def expected_squared_step(self, positions: ArrayLike) -> float:
        x = np.asarray(positions, dtype=float)
        if x.shape != (self.process.size,):
            raise ValueError("positions has the wrong shape")
        squared_displacement = (x[:, None] - x[None, :]) ** 2
        total = 0.0
        for initial, weight in enumerate(self.process.initial_distribution):
            if weight == 0:
                continue
            for coupling in self._couplings(initial):
                total += weight * float(np.sum(coupling * squared_displacement))
        return total


def empirical_marginals(paths: NDArray[np.int64], size: int) -> NDArray[np.float64]:
    """Column-wise empirical distributions for sampled paths."""
    if paths.ndim != 2:
        raise ValueError("paths must be a two-dimensional integer array")
    if len(paths) == 0:
        raise ValueError("cannot estimate marginals from zero paths")
    result = np.zeros((paths.shape[1], size))
    for column in range(paths.shape[1]):
        result[column] = np.bincount(paths[:, column], minlength=size) / len(paths)
    return result
