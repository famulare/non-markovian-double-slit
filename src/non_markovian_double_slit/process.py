"""Positive-law representation of a finite indivisible stochastic process.

This module intentionally has no Hamiltonian or complex-amplitude dependency.
It is the boundary between construction of a law and its classical stochastic
interpretation.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from numpy.typing import ArrayLike, NDArray

ProbabilityVector = NDArray[np.float64]
StochasticMatrix = NDArray[np.float64]


def _probability_vector(values: ArrayLike, size: int, name: str) -> ProbabilityVector:
    if np.iscomplexobj(values):
        raise ValueError(f"{name} must be real")
    vector = np.asarray(values, dtype=float).copy()
    if vector.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},), got {vector.shape}")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} contains non-finite values")
    if np.min(vector) < -1e-13:
        raise ValueError(f"{name} contains a negative probability")
    np.maximum(vector, 0.0, out=vector)
    if not np.isclose(vector.sum(), 1.0, rtol=0, atol=1e-11):
        raise ValueError(f"{name} must sum to one; got {vector.sum():.16g}")
    vector.setflags(write=False)
    return vector


def _stochastic_matrix(values: ArrayLike, size: int, name: str) -> StochasticMatrix:
    if np.iscomplexobj(values):
        raise ValueError(f"{name} must be real")
    matrix = np.asarray(values, dtype=float).copy()
    if matrix.shape != (size, size):
        raise ValueError(f"{name} must have shape ({size}, {size}), got {matrix.shape}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} contains non-finite values")
    if np.min(matrix) < -1e-13:
        raise ValueError(f"{name} contains a negative transition probability")
    np.maximum(matrix, 0.0, out=matrix)
    if not np.allclose(matrix.sum(axis=0), 1.0, rtol=0, atol=1e-11):
        raise ValueError(f"{name} is not column-stochastic")
    matrix.setflags(write=False)
    return matrix


@dataclass(frozen=True)
class FiniteIndivisibleProcess:
    """Barandes's finite process tuple ``(C, T, T0, Gamma, p, A)``.

    ``laws[(target, conditioning)]`` stores the column-stochastic matrix
    ``Gamma(target <- conditioning)``. Conditioning times must belong to
    ``conditioning_times``. In Barandes's tuple, ``A`` is the maximal
    commutative algebra of real random variables ``A: C x T -> R``. It is
    represented here by accepting value vectors when expectations are needed.
    The ordinary event algebra at any fixed time is separately the power set of
    the finite configuration set.

    The object stores no amplitudes. A Hamiltonian-based constructor may create
    these arrays, but downstream evolution and event probabilities are purely
    real and nonnegative.
    """

    configurations: tuple[str, ...]
    times: tuple[float, ...]
    conditioning_times: frozenset[float]
    laws: Mapping[tuple[float, float], ArrayLike]
    initial_time: float
    initial_distribution: ArrayLike

    def __post_init__(self) -> None:
        size = len(self.configurations)
        if size == 0 or len(set(self.configurations)) != size:
            raise ValueError("configurations must be nonempty and unique")
        if len(set(self.times)) != len(self.times):
            raise ValueError("times must be unique")
        time_set = set(self.times)
        if self.initial_time not in time_set:
            raise ValueError("initial_time must belong to times")
        if self.initial_time not in self.conditioning_times:
            raise ValueError("initial_time must be a conditioning time")
        if not self.conditioning_times <= time_set:
            raise ValueError("conditioning_times must be a subset of times")

        checked: dict[tuple[float, float], StochasticMatrix] = {}
        for (target, conditioning), values in self.laws.items():
            if target not in time_set or conditioning not in self.conditioning_times:
                raise ValueError(f"invalid law key {(target, conditioning)}")
            checked[(target, conditioning)] = _stochastic_matrix(
                values, size, f"Gamma({target} <- {conditioning})"
            )

        identity = np.eye(size)
        for conditioning in self.conditioning_times:
            key = (conditioning, conditioning)
            if key not in checked:
                raise ValueError(f"missing trivial law Gamma{key}")
            if not np.allclose(checked[key], identity, rtol=0, atol=1e-11):
                raise ValueError(f"Gamma{key} must be the identity")

        expected_keys = {
            (target, conditioning)
            for target in self.times
            for conditioning in self.conditioning_times
        }
        missing_keys = expected_keys - checked.keys()
        if missing_keys:
            raise ValueError(
                f"laws must define Gamma on T x T0; missing {sorted(missing_keys)}"
            )
        for target in self.times:
            for first in self.conditioning_times:
                for division in self.conditioning_times:
                    direct = checked[(target, first)]
                    divided = checked[(target, division)] @ checked[(division, first)]
                    if not np.allclose(direct, divided, rtol=0, atol=1e-11):
                        raise ValueError(
                            "laws violate composition between conditioning times: "
                            f"Gamma({target} <- {first}) != "
                            f"Gamma({target} <- {division}) "
                            f"Gamma({division} <- {first})"
                        )

        initial = _probability_vector(
            self.initial_distribution, size, "initial_distribution"
        )
        object.__setattr__(self, "laws", MappingProxyType(checked))
        object.__setattr__(self, "initial_distribution", initial)

    @property
    def size(self) -> int:
        return len(self.configurations)

    @property
    def event_algebra_description(self) -> str:
        return f"power set of C ({2**self.size} events)"

    @property
    def random_variable_algebra_description(self) -> str:
        """Description of the ``A`` component of Barandes's tuple."""
        return "maximal commutative algebra of maps A: C x T -> R"

    def law(self, target: float, conditioning: float | None = None) -> StochasticMatrix:
        """Return ``Gamma(target <- conditioning)`` without extrapolation."""
        if conditioning is None:
            conditioning = self.initial_time
        try:
            return self.laws[(target, conditioning)]
        except KeyError as error:
            raise KeyError(
                f"no law conditions at {conditioning} and targets {target}; "
                "indivisible processes do not imply missing conditional laws"
            ) from error

    def marginal(
        self,
        target: float,
        *,
        conditioning: float | None = None,
        distribution: ArrayLike | None = None,
    ) -> ProbabilityVector:
        """Apply a stored positive law to a probability distribution."""
        if distribution is None:
            if conditioning not in (None, self.initial_time):
                raise ValueError(
                    "a noninitial conditioning time requires its own distribution"
                )
            source = self.initial_distribution
        else:
            source = _probability_vector(distribution, self.size, "distribution")
        result = self.law(target, conditioning) @ source
        result.setflags(write=False)
        return result

    def event_probability(
        self,
        target: float,
        event: Iterable[int],
        *,
        conditioning: float | None = None,
        distribution: ArrayLike | None = None,
    ) -> float:
        """Probability of an event in the power-set event algebra."""
        indices = np.unique(np.fromiter(event, dtype=int))
        if np.any((indices < 0) | (indices >= self.size)):
            raise IndexError("event contains an invalid configuration index")
        return float(
            self.marginal(target, conditioning=conditioning, distribution=distribution)[
                indices
            ].sum()
        )

    def expectation(
        self,
        target: float,
        values: ArrayLike,
        *,
        conditioning: float | None = None,
        distribution: ArrayLike | None = None,
    ) -> float:
        """Expectation of a real random variable at one target time."""
        if np.iscomplexobj(values):
            raise ValueError("observable values must be real")
        magnitudes = np.asarray(values, dtype=float)
        if magnitudes.shape != (self.size,) or not np.all(np.isfinite(magnitudes)):
            raise ValueError(f"values must be a finite vector of shape ({self.size},)")
        return float(
            magnitudes
            @ self.marginal(
                target, conditioning=conditioning, distribution=distribution
            )
        )

    def divisibility_residual(
        self,
        target: float,
        intermediate: float,
        candidate_later_law: ArrayLike,
        *,
        conditioning: float | None = None,
    ) -> StochasticMatrix:
        """Return direct law minus a proposed divided composition.

        ``candidate_later_law`` is supplied explicitly because a law conditioned
        at ``intermediate`` need not be part of this process.
        """
        if conditioning is None:
            conditioning = self.initial_time
        later = _stochastic_matrix(
            candidate_later_law, self.size, "candidate_later_law"
        )
        return self.law(target, conditioning) - later @ self.law(
            intermediate, conditioning
        )


def stochastic_error(matrix: ArrayLike) -> float:
    """Maximum negativity or column-normalization error."""
    values = np.asarray(matrix, dtype=float)
    negativity = max(0.0, float(-values.min()))
    normalization = float(np.max(np.abs(values.sum(axis=0) - 1.0)))
    return max(negativity, normalization)
