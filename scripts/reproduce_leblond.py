#!/usr/bin/env python3
"""Reproduce and correct the indexing in LeBlond arXiv:2511.02863v1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_markovian_double_slit.diagnostics import central_fringe_spacing
from non_markovian_double_slit.fresnel import (
    LeBlondResult,
    fraunhofer_double_slit_density,
    leblond_reproduction,
    paraxial_fringe_spacing,
)


def summarize(result: LeBlondResult) -> dict[str, float]:
    effective_separation = float(
        result.upper_points.mean() - result.lower_points.mean()
    )
    expected = paraxial_fringe_spacing(1.23e-10, 1.0, effective_separation)
    observed = central_fringe_spacing(
        result.screen,
        result.coherent_density,
        expected_spacing=expected,
        half_window=0.06,
    )
    midpoint = (result.upper_points.mean() + result.lower_points.mean()) / 2
    analytic = fraunhofer_double_slit_density(
        result.screen - midpoint,
        wavelength=1.23e-10,
        distance=1.0,
        slit_width=1.5e-9,
        slit_separation=effective_separation,
    )
    return {
        "lower_slit_center_m": float(result.lower_points.mean()),
        "upper_slit_center_m": float(result.upper_points.mean()),
        "effective_separation_m": effective_separation,
        "expected_spacing_for_effective_separation_m": expected,
        "observed_spacing_m": observed,
        "fraunhofer_peak_spacing_m": central_fringe_spacing(
            result.screen, analytic, expected_spacing=expected, half_window=0.06
        ),
        "fraunhofer_max_density_error_over_peak": float(
            np.max(np.abs(analytic - result.coherent_density)) / analytic.max()
        ),
        "displayed_probability_mass": result.displayed_mass,
        "plot_coordinate_offset_m": float(
            (result.plotted_screen - result.screen).mean()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", type=int, default=2000)
    parser.add_argument(
        "--output", type=Path, default=Path("outputs/leblond_reproduction.png")
    )
    parser.add_argument(
        "--metrics", type=Path, default=Path("outputs/leblond_metrics.json")
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metrics.parent.mkdir(parents=True, exist_ok=True)

    published = leblond_reproduction(published_upper_indexing=True, size=args.size)
    corrected = leblond_reproduction(published_upper_indexing=False, size=args.size)
    figure, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, sharey=True)
    for axis, result, title in (
        (axes[0], published, "Published Listing 2 indexing"),
        (axes[1], corrected, "Corrected upper-slit local indexing"),
    ):
        axis.plot(result.screen, result.coherent_density, label="coherent")
        axis.plot(result.screen, result.recorded_density, label="which-path mixture")
        axis.set(title=title, ylabel="probability density")
        axis.legend()
    axes[1].set_xlabel("screen position (m)")
    figure.tight_layout()
    figure.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(figure)

    metrics = {"published": summarize(published), "corrected": summarize(corrected)}
    args.metrics.write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
