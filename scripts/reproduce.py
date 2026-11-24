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
    monotone_transport_coupling,
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
    geometry = dict(size=size, length=80.0, slit_separation=8.0, slit_width=2.0)
    model = FiniteDoubleSlitModel.create(**geometry)
    opposite_phase = FiniteDoubleSlitModel.create(**geometry, relative_phase=np.pi)
    elapsed = 4.0
    upper, lower = model.slit_components(elapsed)
    coherent = screen_pattern(upper, lower)
    pi_shifted = opposite_phase.probability(elapsed)
    direct = model.probability(elapsed)
    slit_mass = np.abs(model.slit_state) ** 2
    opposite_slit_mass = np.abs(opposite_phase.slit_state) ** 2
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, first, second, window, title in (
        (axes[0], slit_mass, opposite_slit_mass, 8, "At the slits"),
        (axes[1], direct, pi_shifted, 12, "After the same free flight"),
    ):
        ax.plot(model.grid, first, color="#155773", lw=2, label=r"$\phi=0$")
        ax.plot(model.grid, second, "--", color="#b04427", label=r"$\phi=\pi$")
        ax.set(title=title, xlabel="transverse position", xlim=(-window, window))
        ax.legend(frameon=False)
    axes[0].set_ylabel("probability per cell")
    _save(figure, output, "finite_contexts.png")

    divided = model.divided_law(elapsed, 0.0)[:, model.source_index]
    recorded = retained_record_pattern(upper, lower)
    figure, axis = plt.subplots(figsize=(9, 3.8))
    axis.plot(model.grid, direct, color="#155773", label="no marker")
    axis.plot(model.grid, recorded, "--", color="#b04427", label="opening recorded")
    axis.plot(model.grid, divided, ":", color="#30383a", lw=2, label="cell recorded")
    axis.set(xlabel="transverse position", ylabel="probability per cell")
    axis.legend(frameon=False)
    _save(figure, output, "indivisibility.png")
    residual = direct - divided

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

    # Two source distributions passed through the very same completed P.
    # The second samples a source index uniformly from the aperture support.
    aperture_indices = (model.upper_mode != 0) | (model.lower_mode != 0)
    mixed_source = aperture_indices / aperture_indices.sum()
    early_law = process.law(0.0)
    late_law = process.law(elapsed)
    right = model.grid > 0
    upper_slit = model.upper_mode != 0
    rank_joint = monotone_transport_coupling(slit_mass, direct)
    upper_and_right = rank_joint[np.ix_(right, upper_slit)].sum()
    other_and_right = rank_joint[np.ix_(right, ~upper_slit)].sum()
    right_mass = direct[right].sum()
    upper_mass = slit_mass[upper_slit].sum()
    marked_joint = model.upper_weight * np.abs(upper[right]) ** 2

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
        "restart_tests": {
            "phase_pair_slit_tv": total_variation_distance(
                slit_mass, opposite_slit_mass
            ),
            "phase_pair_screen_tv": total_variation_distance(direct, pi_shifted),
            "center_bin_phase_zero": float(direct[size // 2]),
            "center_bin_phase_pi": float(pi_shifted[size // 2]),
            "fixed_model_slit_tv": total_variation_distance(
                early_law[:, model.source_index], early_law @ mixed_source
            ),
            "fixed_model_screen_tv": total_variation_distance(
                direct, late_law @ mixed_source
            ),
            "coherent_vs_slit_record_tv": total_variation_distance(direct, recorded),
            "slit_vs_cell_record_tv": total_variation_distance(recorded, divided),
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
            "unmarked_right_screen_probability": float(right_mass),
            "upper_given_right_independent": float(upper_mass),
            "upper_given_right_shared_rank": float(
                upper_and_right / (upper_and_right + other_and_right)
            ),
            "upper_given_right_min": float(
                max(0, upper_mass + right_mass - 1) / right_mass
            ),
            "upper_given_right_max": float(min(upper_mass, right_mass) / right_mass),
            "upper_given_right_with_marker": float(
                marked_joint.sum() / recorded[right].sum()
            ),
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
