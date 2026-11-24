"""Supporting two-dimensional split-step double-slit simulation.

The exact finite model in :mod:`non_markovian_double_slit.finite` carries the
formal stochastic construction. This module checks that its transverse physics
is compatible with a finite real barrier and a screen defined by integrated
probability current.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class Lattice2DConfig:
    nx: int = 320
    nz: int = 320
    x_extent: float = 80.0
    z_extent: float = 80.0
    slit_separation: float = 6.0
    slit_width: float = 2.0
    barrier_thickness: float = 1.0
    barrier_height: float = 80.0
    initial_z: float = -15.0
    sigma_x: float = 6.0
    sigma_z: float = 1.5
    wave_number: float = 5.5
    detector_z: float = 14.0
    dt: float = 0.025
    duration: float = 7.2
    mass: float = 1.0
    hbar: float = 1.0

    def validate(self) -> None:
        if self.nx < 32 or self.nz < 32:
            raise ValueError("each lattice dimension must be at least 32")
        positive = (
            self.x_extent,
            self.z_extent,
            self.slit_separation,
            self.slit_width,
            self.barrier_thickness,
            self.barrier_height,
            self.sigma_x,
            self.sigma_z,
            self.wave_number,
            self.dt,
            self.duration,
            self.mass,
            self.hbar,
        )
        if min(positive) <= 0:
            raise ValueError(
                "lengths, scales, duration, mass, and hbar must be positive"
            )
        if self.slit_separation <= self.slit_width:
            raise ValueError("slits must be disjoint")
        if abs(self.initial_z) >= self.z_extent / 2:
            raise ValueError("initial packet center is outside the lattice")
        if abs(self.detector_z) >= self.z_extent / 2:
            raise ValueError("detector is outside the lattice")
        if self.initial_z >= 0 or self.detector_z <= 0:
            raise ValueError("the wall must lie between source and detector")


def config_for_profile(profile: Literal["quick", "reference"]) -> Lattice2DConfig:
    if profile == "quick":
        return Lattice2DConfig(nx=240, nz=320, dt=0.03)
    if profile == "reference":
        return Lattice2DConfig(nx=400, nz=400, dt=0.02)
    raise ValueError(f"unknown profile {profile!r}")


@dataclass(frozen=True)
class FluxDiagnostics:
    net_flux: float
    positive_flux: float
    backward_flux: float
    backward_fraction: float
    negative_bin_mass: float


@dataclass(frozen=True)
class Lattice2DResult:
    config: Lattice2DConfig
    x: NDArray[np.float64]
    z: NDArray[np.float64]
    potential: NDArray[np.float64]
    final_state: NDArray[np.complex128]
    integrated_current: NDArray[np.float64]
    backward_flux: float
    initial_norm: float
    final_norm: float
    max_boundary_probability: float
    snapshots: Mapping[float, NDArray[np.float64]]

    @property
    def norm_error(self) -> float:
        return abs(self.final_norm - self.initial_norm)

    def flux_diagnostics(self) -> FluxDiagnostics:
        dx = float(self.x[1] - self.x[0])
        net = float(self.integrated_current.sum() * dx)
        positive = net + self.backward_flux
        negative_bins = float(np.maximum(-self.integrated_current, 0).sum() * dx)
        return FluxDiagnostics(
            net_flux=net,
            positive_flux=positive,
            backward_flux=self.backward_flux,
            backward_fraction=self.backward_flux / positive if positive > 0 else np.inf,
            negative_bin_mass=negative_bins,
        )

    def screen_probability(
        self,
        *,
        max_backward_fraction: float = 1e-5,
        max_boundary_probability: float = 1e-3,
    ) -> NDArray[np.float64]:
        """Return normalized detector-bin probabilities after explicit checks.

        Tiny negative current bins are set to zero only after their total mass is
        reported by :meth:`flux_diagnostics` and checked against the threshold.
        """
        diagnostics = self.flux_diagnostics()
        if diagnostics.net_flux <= 0:
            raise RuntimeError("detector has nonpositive net integrated flux")
        if diagnostics.backward_fraction > max_backward_fraction:
            raise RuntimeError(
                "backward detector current exceeds threshold: "
                f"{diagnostics.backward_fraction:.3g} > {max_backward_fraction:.3g}"
            )
        if self.max_boundary_probability > max_boundary_probability:
            raise RuntimeError(
                "boundary contamination exceeds threshold: "
                f"{self.max_boundary_probability:.3g} > {max_boundary_probability:.3g}"
            )
        bin_mass = np.maximum(self.integrated_current, 0) * (self.x[1] - self.x[0])
        return bin_mass / bin_mass.sum()


def _coordinates(
    config: Lattice2DConfig,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    dx = config.x_extent / config.nx
    dz = config.z_extent / config.nz
    x = (np.arange(config.nx) - config.nx // 2) * dx
    z = (np.arange(config.nz) - config.nz // 2) * dz
    return x, z


def _barrier(
    config: Lattice2DConfig,
    x: NDArray[np.float64],
    z: NDArray[np.float64],
) -> NDArray[np.float64]:
    in_wall = np.abs(z[:, None]) <= config.barrier_thickness / 2
    in_upper = np.abs(x[None, :] - config.slit_separation / 2) <= (
        config.slit_width / 2
    )
    in_lower = np.abs(x[None, :] + config.slit_separation / 2) <= (
        config.slit_width / 2
    )
    return config.barrier_height * (in_wall & ~(in_upper | in_lower))


def _initial_state(
    config: Lattice2DConfig,
    x: NDArray[np.float64],
    z: NDArray[np.float64],
) -> NDArray[np.complex128]:
    envelope = np.exp(
        -(x[None, :] ** 2) / (4 * config.sigma_x**2)
        - ((z[:, None] - config.initial_z) ** 2) / (4 * config.sigma_z**2)
    )
    state = envelope * np.exp(1j * config.wave_number * z[:, None])
    dx = x[1] - x[0]
    dz = z[1] - z[0]
    return state / np.sqrt(np.sum(np.abs(state) ** 2) * dx * dz)


def _norm(state: NDArray[np.complex128], dx: float, dz: float) -> float:
    return float(np.sum(np.abs(state) ** 2) * dx * dz)


def _boundary_probability(state: NDArray[np.complex128], dx: float, dz: float) -> float:
    margin_z = max(2, state.shape[0] // 32)
    margin_x = max(2, state.shape[1] // 32)
    boundary = np.zeros(state.shape, dtype=bool)
    boundary[:margin_z] = True
    boundary[-margin_z:] = True
    boundary[:, :margin_x] = True
    boundary[:, -margin_x:] = True
    return float(np.sum(np.abs(state[boundary]) ** 2) * dx * dz)


def run_lattice_simulation(config: Lattice2DConfig) -> Lattice2DResult:
    """Run unitary split-step evolution and accumulate signed screen current."""
    config.validate()
    x, z = _coordinates(config)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    potential = _barrier(config, x, z)
    state = _initial_state(config, x, z)
    initial_norm = _norm(state, dx, dz)

    kx = 2 * np.pi * np.fft.fftfreq(config.nx, d=dx)
    kz = 2 * np.pi * np.fft.fftfreq(config.nz, d=dz)
    steps = int(np.ceil(config.duration / config.dt))
    dt = config.duration / steps
    kinetic_phase = np.exp(
        -0.5j * config.hbar * (kz[:, None] ** 2 + kx[None, :] ** 2) * dt / config.mass
    )
    potential_half_phase = np.exp(-0.5j * potential * dt / config.hbar)
    detector_index = int(np.argmin(np.abs(z - config.detector_z)))
    # Row of F† diag(ik_z) F at the detector. A centered finite difference
    # would measure sin(k_z dz)/dz for a plane wave of wave number k_z.
    weights = np.roll(np.fft.fft(1j * kz) / config.nz, detector_index)

    def detector_current() -> NDArray[np.float64]:
        derivative = np.einsum("z,zx->x", weights, state)
        return (config.hbar / config.mass) * np.imag(
            state[detector_index].conj() * derivative
        )

    snapshot_steps = {max(1, round(f * steps)) for f in (0.35, 0.55, 0.8, 1.0)}
    snapshots: dict[float, NDArray[np.float64]] = {}
    # Trapezoidal quadrature, splitting by sign before summing over time.
    initial_current = detector_current()
    integrated_current = 0.5 * dt * initial_current
    backward_flux = 0.5 * dt * dx * np.maximum(-initial_current, 0).sum()
    max_boundary = _boundary_probability(state, dx, dz)

    for step in range(1, steps + 1):
        state *= potential_half_phase
        momentum_state = np.fft.fft2(state, norm="ortho")
        state = np.fft.ifft2(momentum_state * kinetic_phase, norm="ortho")
        state *= potential_half_phase

        current = detector_current()
        weight = 0.5 * dt if step == steps else dt
        integrated_current += weight * current
        backward_flux += weight * dx * np.maximum(-current, 0).sum()
        max_boundary = max(max_boundary, _boundary_probability(state, dx, dz))
        if step in snapshot_steps:
            snapshots[step * dt] = np.abs(state) ** 2

    final_norm = _norm(state, dx, dz)
    return Lattice2DResult(
        config=config,
        x=x,
        z=z,
        potential=potential,
        final_state=state,
        integrated_current=integrated_current,
        backward_flux=float(backward_flux),
        initial_norm=initial_norm,
        final_norm=final_norm,
        max_boundary_probability=max_boundary,
        snapshots=MappingProxyType(snapshots),
    )
