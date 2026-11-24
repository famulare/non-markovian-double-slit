#!/usr/bin/env python3
"""Generate the main finite-process and 2D double-slit demonstrations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from non_markovian_double_slit.diagnostics import (
    central_fringe_spacing,
    normalized_correlation,
    periodic_bin_masses,
    total_variation_distance,
)
from non_markovian_double_slit.finite import (
    FiniteDoubleSlitModel,
    unitarity_error,
)
from non_markovian_double_slit.lattice2d import (
    config_for_profile,
    run_lattice_simulation,
)
from non_markovian_double_slit.process import stochastic_error
from non_markovian_double_slit.realizers import (
    IndependentPathRealizer,
    TransportPathRealizer,
)
from non_markovian_double_slit.which_path import (
    coherent_uncomputation_pattern,
    local_reset_pattern,
    quantum_eraser_patterns,
    retained_record_pattern,
    screen_pattern,
)


def _save(figure: plt.Figure, output: Path, name: str) -> None:
    figure.tight_layout()
    figure.savefig(output / name, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _barandes_two_state() -> dict[str, list[float]]:
    hadamard = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2)
    source = np.array([1.0, 0.0])
    direct = np.abs(hadamard @ hadamard) ** 2 @ source
    divided = (np.abs(hadamard) ** 2 @ np.abs(hadamard) ** 2) @ source
    return {
        "direct": direct.tolist(),
        "divided": divided.tolist(),
        "interference_residual": (direct - divided).tolist(),
    }


def finite_demonstrations(
    profile: str, output: Path
) -> tuple[FiniteDoubleSlitModel, dict[str, object]]:
    size = 201 if profile == "quick" else 401
    model = FiniteDoubleSlitModel.create(
        size=size,
        length=80.0,
        slit_separation=8.0,
        slit_width=2.0,
    )
    elapsed = 4.0
    upper, lower = model.slit_components(elapsed)
    coherent = screen_pattern(upper, lower)
    upper_only = np.abs(upper) ** 2
    lower_only = np.abs(lower) ** 2
    pi_shifted = screen_pattern(upper, lower, relative_phase=np.pi)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(model.grid, coherent, label="both slits")
    axes[0].plot(model.grid, 0.5 * upper_only, "--", label="upper contribution")
    axes[0].plot(model.grid, 0.5 * lower_only, "--", label="lower contribution")
    axes[0].set(
        title="Slit contexts",
        xlabel="screen configuration x",
        ylabel="probability",
    )
    axes[0].legend()
    axes[1].plot(model.grid, coherent, label="phase 0")
    axes[1].plot(model.grid, pi_shifted, label="phase pi")
    axes[1].set(title="Relative phase moves the fringes", xlabel="x")
    axes[1].legend()
    _save(figure, output, "finite_contexts.png")

    direct = model.probability(elapsed)
    divided = model.divided_law(elapsed, 0.0)[:, model.source_index]
    residual = direct - divided
    figure, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    axes[0].plot(model.grid, direct, label="direct Gamma(t <- source)")
    axes[0].plot(model.grid, divided, label="position recorded to one cell")
    axes[0].set(
        ylabel="probability",
        title="Recording the aperture cell destroys within-slit coherence too",
    )
    axes[0].legend()
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].plot(model.grid, residual, color="tab:red")
    axes[1].set(xlabel="screen configuration x", ylabel="direct - divided")
    _save(figure, output, "indivisibility.png")

    recorded = retained_record_pattern(upper, lower)
    local_reset = local_reset_pattern(upper, lower)
    uncomputed = coherent_uncomputation_pattern(upper, lower)
    partial = screen_pattern(upper, lower, record_overlap=0.5)
    figure, axis = plt.subplots(figsize=(9, 4))
    axis.plot(model.grid, coherent, label="no record / coherent uncomputation")
    axis.plot(model.grid, partial, label="record overlap 0.5")
    axis.plot(model.grid, recorded, label="record retained / local reset")
    axis.set(
        title="Which-path information controls the cross term",
        xlabel="screen configuration x",
        ylabel="probability",
    )
    axis.legend()
    _save(figure, output, "which_path.png")

    eraser = quantum_eraser_patterns(upper, lower)
    figure, axis = plt.subplots(figsize=(9, 4))
    axis.plot(model.grid, eraser.plus / eraser.plus.sum(), label="conditional +")
    axis.plot(model.grid, eraser.minus / eraser.minus.sum(), label="conditional -")
    axis.plot(model.grid, eraser.aggregate, color="black", lw=2, label="unconditioned")
    axis.set(
        title="Quantum eraser: conditional fringe and antifringe",
        xlabel="screen configuration x",
        ylabel="conditional or aggregate probability",
    )
    axis.legend()
    _save(figure, output, "quantum_eraser.png")

    rng = np.random.default_rng(20260721)
    figure, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharey=True)
    for axis, count in zip(axes, (50, 500, 5000), strict=True):
        events = rng.choice(model.size, size=count, p=direct)
        empirical = np.bincount(events, minlength=model.size) / count
        axis.step(model.grid, empirical, where="mid", label="events")
        axis.plot(model.grid, direct, color="black", lw=1.3, label="Gamma column")
        axis.set(title=f"{count:,} detections", xlabel="x")
    axes[0].set_ylabel("relative frequency")
    axes[-1].legend()
    _save(figure, output, "event_buildup.png")

    target_times = (0.0, 1.0, 2.0, elapsed)
    process = model.process(target_times)
    independent = IndependentPathRealizer(process, target_times)
    transport = TransportPathRealizer(process, target_times)
    independent_paths = independent.sample(35, rng=np.random.default_rng(7))
    transport_paths = transport.sample(35, rng=np.random.default_rng(7))
    display_times = np.array([model.source_time, *target_times])
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for axis, paths, title in (
        (axes[0], independent_paths, "conditionally independent"),
        (axes[1], transport_paths, "monotone transport"),
    ):
        for path in paths:
            axis.plot(display_times, model.grid[path], alpha=0.28, lw=0.8)
        axis.set(title=title, xlabel="target time")
    axes[0].set_ylabel("configuration x")
    _save(figure, output, "path_realizers.png")

    metrics: dict[str, object] = {
        "barandes_two_state": _barandes_two_state(),
        "finite": {
            "size": model.size,
            "preparation_unitarity_error": unitarity_error(model.preparation),
            "evolution_unitarity_error": unitarity_error(model.evolution(elapsed)),
            "gamma_stochastic_error": stochastic_error(model.gamma(elapsed)),
            "direct_divided_total_variation": total_variation_distance(direct, divided),
            "max_abs_interference_residual": float(np.max(np.abs(residual))),
        },
        "which_path": {
            "local_reset_vs_retained_max_error": float(
                np.max(np.abs(local_reset - recorded))
            ),
            "uncomputation_vs_coherent_max_error": float(
                np.max(np.abs(uncomputed - coherent))
            ),
            "eraser_aggregate_vs_recorded_max_error": float(
                np.max(np.abs(eraser.aggregate - recorded))
            ),
        },
        "path_realizers": {
            "independent_expected_squared_step": independent.expected_squared_step(
                model.grid
            ),
            "transport_expected_squared_step": transport.expected_squared_step(
                model.grid
            ),
        },
    }
    return model, metrics


def lattice_demonstration(profile: str, output: Path) -> dict[str, object]:
    config = config_for_profile(profile)
    result = run_lattice_simulation(config)
    screen = result.screen_probability()
    flux = result.flux_diagnostics()
    expected_spacing = (
        2 * np.pi * config.detector_z / (config.wave_number * config.slit_separation)
    )
    observed_spacing = central_fringe_spacing(
        result.x, screen, expected_spacing=expected_spacing
    )

    snapshot_time = sorted(result.snapshots)[1]
    snapshot = result.snapshots[snapshot_time]
    figure, axes = plt.subplots(1, 2, figsize=(12, 4))
    image = axes[0].imshow(
        np.log10(snapshot + 1e-12),
        origin="lower",
        aspect="auto",
        extent=(result.x[0], result.x[-1], result.z[0], result.z[-1]),
        cmap="magma",
        vmin=-9,
    )
    axes[0].contour(
        result.x,
        result.z,
        result.potential,
        levels=[config.barrier_height / 2],
        colors=["cyan"],
        linewidths=0.8,
    )
    axes[0].axhline(config.detector_z, color="white", ls="--", lw=0.8)
    axes[0].set(
        title=f"2D probability density at t={snapshot_time:.2f}",
        xlabel="transverse x",
        ylabel="longitudinal z",
    )
    axes[0].legend(
        handles=[
            Line2D([0], [0], color="cyan", lw=0.8, label="finite slit wall"),
            Line2D(
                [0],
                [0],
                color="white",
                ls="--",
                lw=0.8,
                label="detector",
            ),
        ],
        loc="upper right",
    )
    figure.colorbar(image, ax=axes[0], label="log10 probability density")
    axes[1].plot(result.x, screen)
    axes[1].set(
        title="Normalized time-integrated detector current",
        xlabel="transverse x",
        ylabel="detector-bin probability",
        xlim=(-15, 15),
    )
    _save(figure, output, "lattice2d.png")

    transverse = FiniteDoubleSlitModel.create(
        size=401 if profile == "quick" else 601,
        length=config.x_extent,
        slit_separation=config.slit_separation,
        slit_width=config.slit_width,
        mass=config.mass,
        hbar=config.hbar,
    )
    velocity_z = config.hbar * config.wave_number / config.mass
    flight_time = config.detector_z / velocity_z
    exact = periodic_bin_masses(
        transverse.grid, transverse.probability(flight_time), result.x
    )
    exact_spacing = central_fringe_spacing(
        result.x, exact, expected_spacing=expected_spacing
    )
    figure, axis = plt.subplots(figsize=(9, 4))
    axis.plot(result.x, screen, label="2D integrated current")
    axis.plot(result.x, exact, "--", label="finite transverse law")
    axis.set(
        title="Paraxial cross-check",
        xlabel="transverse x",
        ylabel="normalized detector-bin probability",
        xlim=(-15, 15),
    )
    axis.legend()
    _save(figure, output, "exact_vs_2d.png")

    return {
        "shape": list(result.final_state.shape),
        "norm_error": result.norm_error,
        "max_boundary_probability": result.max_boundary_probability,
        "net_detector_flux": flux.net_flux,
        "backward_flux_fraction": flux.backward_fraction,
        "negative_integrated_bin_mass": flux.negative_bin_mass,
        "paraxial_expected_spacing": expected_spacing,
        "observed_spacing": observed_spacing,
        "exact_finite_spacing": exact_spacing,
        "exact_2d_correlation": normalized_correlation(screen, exact),
        "exact_2d_total_variation": total_variation_distance(screen, exact),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("quick", "reference"), default="quick")
    parser.add_argument("--output", type=Path, default=Path("outputs/figures"))
    parser.add_argument("--metrics", type=Path, default=Path("outputs/metrics.json"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.metrics.parent.mkdir(parents=True, exist_ok=True)

    _, metrics = finite_demonstrations(args.profile, args.output)
    metrics["lattice2d"] = lattice_demonstration(args.profile, args.output)
    args.metrics.write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
