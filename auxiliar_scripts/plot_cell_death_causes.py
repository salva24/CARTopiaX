"""
Copyright 2026 compiler-research.org, Salvador de la Torre Gonzalez

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
    SPDX-License-Identifier: Apache-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This file plots, in a single graph and at a given minute, the radial profile of
the tumor cell deaths for each cause (lack of oxygen, lack of glucose, random
natural causes and CAR-T kill) and prints the total number of deaths of each
cause (summed over all radii), using output/data_dependent_on_radius_tumor.csv. The CSV has one
column per attribute and ring, named "<attribute>_<r_in>_to_<r_out>".
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analysis_radius_attribute import load_profile

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
CSV_PATH = BASE_DIR / "output" / "data_dependent_on_radius_tumor.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "out"  # destination folder

# Simulation minute to plot (must exist in the "total_minutes" column)
MINUTE = 0

# Death causes: column prefix -> legend label
DEATH_CAUSES = {
    "dead_tumor_cells_lack_of_oxygen_radius": "Lack of oxygen",
    "dead_tumor_cells_lack_of_glucose_radius": "Lack of glucose",
    "dead_tumor_cells_random_natural_causes_radius": "Random natural causes",
    "dead_tumor_cells_cartcell_kill_radius": "CAR-T kill",
}

# True: divide counts by the ring area, so the value is expressed per um^2
# instead of as the total amount inside the ring
APPLY_RADIAL_NORMALIZATION = True

# Y axis limits (None = automatic)
Y_MIN = 0
Y_MAX = None


def main():
    fig, ax = plt.subplots(figsize=(9, 5.5))
    radii = None
    total_deaths = {}

    for attribute, label in DEATH_CAUSES.items():
        radii, r_in, r_out, values = load_profile(CSV_PATH, attribute, MINUTE)
        # Total over all rings (raw counts, independent of radius/normalization)
        total_deaths[label] = float(np.sum(values))
        if APPLY_RADIAL_NORMALIZATION:
            values = values / (np.pi * (r_out**2 - r_in**2))
        ax.plot(radii, values, marker="o", label=label)

    print(f"Minute: {MINUTE}")
    print("Total deaths by cause (all radii):")
    for label, count in total_deaths.items():
        print(f"  {label}: {count:g}")
    print(f"Radius range with data: {radii[0]:g} to {radii[-1]:g}")
    print(f"Number of points: {len(radii)}")
    print(f"Radial normalization applied: {APPLY_RADIAL_NORMALIZATION}")

    ax.set_xlabel("Radius (um)")
    ax.set_ylabel(
        "Dead tumor cells" + (" (per um^2)" if APPLY_RADIAL_NORMALIZATION else "")
    )
    ax.set_title(f"Tumor cell deaths by cause vs radius - minute {MINUTE}")
    ax.grid(True, alpha=0.3)
    ax.legend()
    if Y_MIN is not None or Y_MAX is not None:
        ax.set_ylim(bottom=Y_MIN, top=Y_MAX)
    fig.tight_layout()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / f"cell_death_causes_vs_radius_minute_{MINUTE}.png"
    fig.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
