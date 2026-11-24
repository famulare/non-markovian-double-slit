import numpy as np

from non_markovian_double_slit.finite import (
    FiniteDoubleSlitModel,
    unitarity_error,
)
from non_markovian_double_slit.process import stochastic_error


def model(phase: float = 0.0) -> FiniteDoubleSlitModel:
    return FiniteDoubleSlitModel.create(
        size=51,
        length=51.0,
        slit_separation=8.0,
        slit_width=3.0,
        relative_phase=phase,
    )


def test_full_finite_construction_is_unitary_and_unistochastic() -> None:
    double_slit = model()
    assert unitarity_error(double_slit.preparation) < 1e-12
    assert unitarity_error(double_slit.evolution(2.0)) < 1e-12
    assert stochastic_error(double_slit.gamma(2.0)) < 1e-12
    assert np.allclose(
        double_slit.preparation @ double_slit.source_state,
        double_slit.slit_state,
    )
    assert np.isclose(double_slit.probability(2.0).sum(), 1.0)


def test_direct_law_differs_from_divided_stochastic_composition() -> None:
    double_slit = model()
    residual = double_slit.interference_residual(2.0, 0.0)
    assert np.max(np.abs(residual)) > 1e-3
    assert np.min(residual) < 0 < np.max(residual)
    assert np.allclose(residual.sum(axis=0), 0.0, atol=1e-12)


def test_positive_process_contains_only_anchor_conditioned_laws() -> None:
    double_slit = model()
    process = double_slit.process([0.0, 1.0, 2.0])
    assert process.conditioning_times == frozenset({double_slit.source_time})
    assert np.allclose(process.law(2.0), double_slit.gamma(2.0))
    assert np.allclose(process.marginal(2.0), double_slit.probability(2.0))


def test_relative_phase_changes_the_screen_law() -> None:
    zero_phase = model(phase=0.0).probability(2.0)
    pi_phase = model(phase=np.pi).probability(2.0)
    assert not np.allclose(zero_phase, pi_phase)
    assert np.isclose(zero_phase.sum(), 1.0)
    assert np.isclose(pi_phase.sum(), 1.0)


def test_nearly_coincident_source_and_target() -> None:
    from non_markovian_double_slit.finite import unitary_mapping

    a = np.array([1.0, 0.0])
    for angle in (1e-7, 1e-5):
        b = np.array([np.cos(angle), 1j * np.sin(angle)])
        mapping = unitary_mapping(a, b)
        np.testing.assert_allclose(mapping @ a, b, rtol=0, atol=1e-14)
        assert unitarity_error(mapping) < 1e-14


def test_fixed_preparation_has_no_stochastic_restart_at_the_slit_plane() -> None:
    from non_markovian_double_slit.diagnostics import total_variation_distance

    m = FiniteDoubleSlitModel.create(
        size=401, length=80, slit_separation=8, slit_width=2
    )
    process = m.process([0, 4])
    slit, screen = process.law(0), process.law(4)
    cells = (m.upper_mode != 0) | (m.lower_mode != 0)
    mixed = cells / cells.sum()
    distances = [
        total_variation_distance(law[:, m.source_index], law @ mixed)
        for law in (slit, screen)
    ]
    assert abs(distances[0] - 1 / cells.sum()) < 1e-13
    assert distances[1] > distances[0] + 0.68
    # Independent check on all columns. In this completion the earlier law is
    # invertible; its unique later factor violates nonnegativity.
    candidate = np.linalg.solve(slit.T, screen.T).T
    np.testing.assert_allclose(candidate @ slit, screen, atol=1e-13)
    assert candidate.min() < -0.4


def test_interference_residual_can_coexist_with_a_positive_factor() -> None:
    angle = np.pi / 6
    c, s = np.cos(angle), np.sin(angle)
    rotation = np.array([[c, -s], [s, c]])
    first = rotation**2
    final = (rotation @ rotation) ** 2
    assert not np.allclose(first @ first, final)
    np.testing.assert_allclose(first[::-1, :], final)
