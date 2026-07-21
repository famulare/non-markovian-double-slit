import numpy as np

from non_markovian_double_slit.diagnostics import (
    central_fringe_spacing,
    normalized_correlation,
    total_variation_distance,
)
from non_markovian_double_slit.finite import FiniteDoubleSlitModel
from non_markovian_double_slit.lattice2d import (
    config_for_profile,
    run_lattice_simulation,
)


def test_quick_lattice_run_preserves_norm_and_resolves_paraxial_fringes() -> None:
    config = config_for_profile("quick")
    result = run_lattice_simulation(config)
    diagnostics = result.flux_diagnostics()
    screen = result.screen_probability()
    expected = (
        2 * np.pi * config.detector_z / (config.wave_number * config.slit_separation)
    )
    observed = central_fringe_spacing(result.x, screen, expected_spacing=expected)

    assert result.norm_error < 1e-10
    assert result.max_boundary_probability < 1e-3
    assert diagnostics.backward_fraction < 1e-5
    assert np.isclose(screen.sum(), 1.0)
    assert abs(observed / expected - 1) < 0.15

    transverse = FiniteDoubleSlitModel.create(
        size=401,
        length=config.x_extent,
        slit_separation=config.slit_separation,
        slit_width=config.slit_width,
    )
    flight_time = config.detector_z / config.wave_number
    exact = np.interp(
        result.x,
        transverse.grid,
        transverse.probability(flight_time),
        left=0.0,
        right=0.0,
    )
    exact /= exact.sum()
    assert normalized_correlation(screen, exact) > 0.95
    assert total_variation_distance(screen, exact) < 0.15
