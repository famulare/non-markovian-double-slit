"""Continuum Fresnel benchmarks, including LeBlond's published discretization."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

PLANCK = 6.62607015e-34
HBAR = PLANCK / (2 * np.pi)
ELECTRON_MASS = 9.1093837139e-31


def free_kernel(
    screen: ArrayLike,
    aperture: ArrayLike,
    elapsed: float,
    *,
    mass: float = ELECTRON_MASS,
    hbar: float = HBAR,
) -> NDArray[np.complex128]:
    """One-dimensional free-particle kernel ``K(x, t; x', 0)``."""
    if elapsed <= 0 or mass <= 0 or hbar <= 0:
        raise ValueError("elapsed, mass, and hbar must be positive")
    x = np.asarray(screen, dtype=float)[:, None]
    aperture_x = np.asarray(aperture, dtype=float)[None, :]
    prefactor = np.sqrt(mass / (2j * np.pi * hbar * elapsed))
    return prefactor * np.exp(1j * mass * (x - aperture_x) ** 2 / (2 * hbar * elapsed))


def integrate_uniform_slit(
    screen: ArrayLike,
    aperture_points: ArrayLike,
    cell_width: float,
    total_open_width: float,
    elapsed: float,
    *,
    chunk_size: int = 256,
) -> NDArray[np.complex128]:
    """Midpoint quadrature of a uniform aperture amplitude density."""
    x = np.asarray(screen, dtype=float)
    points = np.asarray(aperture_points, dtype=float)
    result = np.empty(x.size, dtype=complex)
    aperture_density = 1 / np.sqrt(total_open_width)
    for start in range(0, x.size, chunk_size):
        stop = min(start + chunk_size, x.size)
        result[start:stop] = (
            free_kernel(x[start:stop], points, elapsed).sum(axis=1)
            * aperture_density
            * cell_width
        )
    return result


@dataclass(frozen=True)
class LeBlondResult:
    screen: NDArray[np.float64]
    plotted_screen: NDArray[np.float64]
    lower_points: NDArray[np.float64]
    upper_points: NDArray[np.float64]
    lower_amplitude: NDArray[np.complex128]
    upper_amplitude: NDArray[np.complex128]
    coherent_density: NDArray[np.float64]
    recorded_density: NDArray[np.float64]
    displayed_mass: float


def leblond_reproduction(
    *,
    published_upper_indexing: bool,
    size: int = 2000,
) -> LeBlondResult:
    """Reproduce LeBlond's v1 quadrature, with optional indexing correction."""
    if size < 4 or size % 2:
        raise ValueError("size must be a positive even integer")
    wavelength = 1.23e-10
    slit_width = 1.5e-9
    slit_separation = 6.15e-9
    distance = 1.0
    velocity = PLANCK / (wavelength * ELECTRON_MASS)
    elapsed = distance / velocity
    screen_min, screen_max = -0.15, 0.15
    screen_width = (screen_max - screen_min) / size
    screen = screen_min + (np.arange(size) + 0.5) * screen_width
    # The plotting tibble in the paper adds one extra screen interval.
    plotted_screen = screen + screen_width

    cells_per_slit = size // 2
    slit_cell_width = slit_width / cells_per_slit
    local = np.arange(cells_per_slit)
    lower_min = -(slit_separation + slit_width) / 2
    upper_min = (slit_separation - slit_width) / 2
    lower_points = lower_min + (local + 0.5) * slit_cell_width
    if published_upper_indexing:
        # Listing 2 uses the global R index ip=N/2+1,...,N in this expression.
        upper_index = local + cells_per_slit
    else:
        upper_index = local
    upper_points = upper_min + (upper_index + 0.5) * slit_cell_width

    lower = integrate_uniform_slit(
        screen, lower_points, slit_cell_width, 2 * slit_width, elapsed
    )
    upper = integrate_uniform_slit(
        screen, upper_points, slit_cell_width, 2 * slit_width, elapsed
    )
    coherent = np.abs(lower + upper) ** 2
    recorded = np.abs(lower) ** 2 + np.abs(upper) ** 2
    return LeBlondResult(
        screen=screen,
        plotted_screen=plotted_screen,
        lower_points=lower_points,
        upper_points=upper_points,
        lower_amplitude=lower,
        upper_amplitude=upper,
        coherent_density=coherent,
        recorded_density=recorded,
        displayed_mass=float(coherent.sum() * screen_width),
    )


def paraxial_fringe_spacing(
    wavelength: float, distance: float, slit_separation: float
) -> float:
    """Small-angle center-to-center fringe spacing ``lambda L / d``."""
    if min(wavelength, distance, slit_separation) <= 0:
        raise ValueError("all arguments must be positive")
    return wavelength * distance / slit_separation
