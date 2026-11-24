import numpy as np

from non_markovian_double_slit.finite import FiniteDoubleSlitModel
from non_markovian_double_slit.which_path import (
    coherent_uncomputation_pattern,
    local_reset_pattern,
    path_distinguishability,
    quantum_eraser_patterns,
    retained_record_pattern,
    screen_pattern,
)


def amplitudes() -> tuple[np.ndarray, np.ndarray]:
    model = FiniteDoubleSlitModel.create(
        size=51, length=51.0, slit_separation=8.0, slit_width=3.0
    )
    return model.slit_components(2.0)


def test_record_overlap_interpolates_between_mixture_and_coherence() -> None:
    upper, lower = amplitudes()
    coherent = screen_pattern(upper, lower, record_overlap=1.0)
    mixture = retained_record_pattern(upper, lower)
    halfway = screen_pattern(upper, lower, record_overlap=0.5)
    assert np.allclose(halfway, 0.5 * (coherent + mixture))
    assert np.isclose(coherent.sum(), 1.0)
    assert np.isclose(mixture.sum(), 1.0)


def test_local_reset_does_not_restore_electron_coherence() -> None:
    upper, lower = amplitudes()
    assert np.allclose(
        local_reset_pattern(upper, lower),
        retained_record_pattern(upper, lower),
    )
    assert not np.allclose(
        local_reset_pattern(upper, lower),
        coherent_uncomputation_pattern(upper, lower),
    )


def test_quantum_eraser_conditionals_sum_to_unconditional_mixture() -> None:
    upper, lower = amplitudes()
    eraser = quantum_eraser_patterns(upper, lower)
    mixture = retained_record_pattern(upper, lower)
    assert np.allclose(eraser.aggregate, mixture)
    assert np.isclose(eraser.plus.sum(), 0.5)
    assert np.isclose(eraser.minus.sum(), 0.5)
    assert not np.allclose(eraser.plus / 0.5, eraser.minus / 0.5)


def test_pure_record_visibility_distinguishability_relation() -> None:
    for overlap in (0.0, 0.25, 0.8, 1.0):
        distinguishability = path_distinguishability(overlap)
        assert np.isclose(distinguishability**2 + overlap**2, 1.0)


def test_distinguishability_includes_unequal_prior_odds() -> None:
    q = 0.8
    for mu in (1.0, 0.0, 0.3j, -0.5):
        upper = np.array([1, 0])
        lower = np.array([mu, np.sqrt(1 - abs(mu) ** 2)])
        difference = q * np.outer(upper, upper) - (1 - q) * np.outer(
            lower, lower.conj()
        )
        assert np.isclose(
            path_distinguishability(mu, upper_weight=q),
            np.abs(np.linalg.eigvalsh(difference)).sum(),
        )


def test_complex_marker_overlap_matches_joint_born_probabilities() -> None:
    upper, lower = amplitudes()
    mu = 0.5j
    q, phase = 0.7, 0.3
    markers = np.array([[1, 0], [mu, np.sqrt(1 - abs(mu) ** 2)]])
    branches = np.column_stack(
        (np.sqrt(q) * upper, np.exp(1j * phase) * np.sqrt(1 - q) * lower)
    )
    joint = branches @ markers
    np.testing.assert_allclose(
        screen_pattern(
            upper, lower, upper_weight=q, relative_phase=phase, record_overlap=mu
        ),
        (np.abs(joint) ** 2).sum(axis=1),
        atol=1e-15,
    )
