"""Positive path-space realizers for prescribed one-time process marginals."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .process import FiniteIndivisibleProcess, _probability_vector


def monotone_transport_coupling(
    source: ArrayLike, target: ArrayLike
) -> NDArray[np.float64]:
    """One-dimensional quadratic-cost optimal coupling.

    Rows index the target state and columns the source state. For ordered
    one-dimensional support, the monotone coupling minimizes every convex
    function of displacement, including squared displacement.
    """
    size = np.size(source)
    source_values = _probability_vector(source, size, "source")
    target_values = _probability_vector(target, size, "target")

    source_edges = np.r_[0.0, np.cumsum(source_values / source_values.sum())]
    target_edges = np.r_[0.0, np.cumsum(target_values / target_values.sum())]
    left = np.maximum(target_edges[:-1, None], source_edges[None, :-1])
    right = np.minimum(target_edges[1:, None], source_edges[None, 1:])
    return np.maximum(right - left, 0.0)


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
        if len(set(self.target_times)) != len(self.target_times):
            raise ValueError("target times must be distinct")
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
        if len(set(self.target_times)) != len(self.target_times):
            raise ValueError("target times must be distinct")
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
