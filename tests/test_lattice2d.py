import numpy as np

from non_markovian_double_slit.diagnostics import (
    central_fringe_spacing,
    normalized_correlation,
    periodic_bin_masses,
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
    assert 0 < diagnostics.backward_fraction < 1e-5
    assert diagnostics.negative_bin_mass < diagnostics.backward_flux
    # Continuity on the effectively open downstream half-domain. The reference
    # plane bisects a grid cell and the far boundary has negligible population.
    detector = np.argmin(np.abs(result.z - config.detector_z))
    density = np.abs(result.final_state) ** 2
    downstream_mass = density[detector + 1 :].sum() + 0.5 * density[detector].sum()
    downstream_mass *= (result.x[1] - result.x[0]) * (result.z[1] - result.z[0])
    assert abs(diagnostics.net_flux - downstream_mass) < 1e-3
    assert np.isclose(screen.sum(), 1.0)
    assert abs(observed / expected - 1) < 0.15

    transverse = FiniteDoubleSlitModel.create(
        size=401,
        length=config.x_extent,
        slit_separation=config.slit_separation,
        slit_width=config.slit_width,
    )
    flight_time = config.detector_z / config.wave_number
    exact = periodic_bin_masses(
        transverse.grid, transverse.probability(flight_time), result.x
    )
    assert normalized_correlation(screen, exact) > 0.95
    assert total_variation_distance(screen, exact) < 0.15


def test_rebin_splits_a_cell_across_the_periodic_seam() -> None:
    actual = periodic_bin_masses([0.5, 1.5, 2.5, 3.5], [0, 0, 0, 1], [0, 2])
    np.testing.assert_allclose(actual, [1, 0])
    actual = periodic_bin_masses([0.5, 1.5, 2.5, 3.5], [0, 0, 0, 1], [0.5, 2.5])
    np.testing.assert_allclose(actual, [0.5, 0.5])
