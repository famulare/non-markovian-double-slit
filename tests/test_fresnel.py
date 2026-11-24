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


def test_far_field_box_slits_predict_the_fresnel_density_in_physical_units() -> None:
    from non_markovian_double_slit.fresnel import fraunhofer_double_slit_density

    for published in (False, True):
        result = leblond_reproduction(published_upper_indexing=published)
        lower, upper = result.lower_points.mean(), result.upper_points.mean()
        analytic = fraunhofer_double_slit_density(
            result.screen - (lower + upper) / 2,
            wavelength=1.23e-10,
            distance=1.0,
            slit_width=1.5e-9,
            slit_separation=upper - lower,
        )
        assert (
            np.max(np.abs(analytic - result.coherent_density)) / analytic.max() < 4e-7
        )
