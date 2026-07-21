import numpy as np
import pytest

from non_markovian_double_slit.process import FiniteIndivisibleProcess


def simple_process() -> FiniteIndivisibleProcess:
    return FiniteIndivisibleProcess(
        configurations=("left", "right"),
        times=(0.0, 1.0),
        conditioning_times=frozenset({0.0}),
        laws={
            (0.0, 0.0): np.eye(2),
            (1.0, 0.0): np.array([[0.7, 0.2], [0.3, 0.8]]),
        },
        initial_time=0.0,
        initial_distribution=np.array([0.25, 0.75]),
    )


def test_process_is_a_positive_law_tuple() -> None:
    process = simple_process()
    assert np.allclose(process.marginal(1.0), [0.325, 0.675])
    assert process.event_probability(1.0, [1]) == pytest.approx(0.675)
    assert process.event_probability(1.0, [1, 1]) == pytest.approx(0.675)
    assert process.expectation(1.0, [-1.0, 1.0]) == pytest.approx(0.35)
    assert process.event_algebra_description == "power set of C (4 events)"
    assert "A: C x T -> R" in process.random_variable_algebra_description
    assert process.law(1.0).dtype == np.dtype(float)


def test_missing_conditioning_law_is_not_inferred() -> None:
    process = simple_process()
    with pytest.raises(KeyError, match="do not imply"):
        process.law(1.0, conditioning=1.0)


def test_process_rejects_nontrivial_law_at_conditioning_time() -> None:
    with pytest.raises(ValueError, match="identity"):
        FiniteIndivisibleProcess(
            configurations=("a", "b"),
            times=(0.0,),
            conditioning_times=frozenset({0.0}),
            laws={(0.0, 0.0): np.full((2, 2), 0.5)},
            initial_time=0.0,
            initial_distribution=[1.0, 0.0],
        )


def test_process_rejects_negative_law() -> None:
    with pytest.raises(ValueError, match="negative"):
        FiniteIndivisibleProcess(
            configurations=("a", "b"),
            times=(0.0, 1.0),
            conditioning_times=frozenset({0.0}),
            laws={
                (0.0, 0.0): np.eye(2),
                (1.0, 0.0): np.array([[1.1, 0.0], [-0.1, 1.0]]),
            },
            initial_time=0.0,
            initial_distribution=[1.0, 0.0],
        )


def test_process_requires_law_for_every_target_conditioning_pair() -> None:
    with pytest.raises(ValueError, match="T x T0"):
        FiniteIndivisibleProcess(
            configurations=("a", "b"),
            times=(0.0, 1.0),
            conditioning_times=frozenset({0.0}),
            laws={(0.0, 0.0): np.eye(2)},
            initial_time=0.0,
            initial_distribution=[1.0, 0.0],
        )


def test_process_enforces_composition_between_division_events() -> None:
    mixing = np.full((2, 2), 0.5)
    with pytest.raises(ValueError, match="violate composition"):
        FiniteIndivisibleProcess(
            configurations=("a", "b"),
            times=(0.0, 1.0),
            conditioning_times=frozenset({0.0, 1.0}),
            laws={
                (0.0, 0.0): np.eye(2),
                (1.0, 1.0): np.eye(2),
                (1.0, 0.0): mixing,
                (0.0, 1.0): mixing,
            },
            initial_time=0.0,
            initial_distribution=[1.0, 0.0],
        )
