import numpy as np

from non_markovian_double_slit.fresnel import (
    leblond_reproduction,
    paraxial_fringe_spacing,
)


def test_leblond_published_upper_index_is_shifted_by_one_slit_width() -> None:
    published = leblond_reproduction(published_upper_indexing=True, size=200)
    corrected = leblond_reproduction(published_upper_indexing=False, size=200)
    shift = published.upper_points.mean() - corrected.upper_points.mean()
    assert np.isclose(shift, 1.5e-9, rtol=1e-12, atol=1e-22)
    assert np.isclose(corrected.upper_points.mean(), 3.075e-9, rtol=1e-12, atol=1e-22)
    assert np.isclose(corrected.lower_points.mean(), -3.075e-9, rtol=1e-12, atol=1e-22)
    assert np.isclose(published.upper_points.mean(), 4.575e-9, rtol=1e-12, atol=1e-22)


def test_leblond_plot_coordinate_is_shifted_by_one_screen_bin() -> None:
    result = leblond_reproduction(published_upper_indexing=True, size=200)
    assert np.allclose(result.plotted_screen - result.screen, 0.3 / 200)


def test_paraxial_electron_fringe_spacing() -> None:
    assert np.isclose(paraxial_fringe_spacing(1.23e-10, 1.0, 6.15e-9), 0.02)
