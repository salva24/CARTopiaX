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

This file plots the radial profile of one attribute of the simulation at a
given minute, using output/data_dependent_on_radius_tumor.csv. The CSV has one
column per attribute and ring, named "<attribute>_<r_in>_to_<r_out>".
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
CSV_PATH = BASE_DIR / "output" / "data_dependent_on_radius_tumor.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "out"  # destination folder

# Simulation minute to plot (must exist in the "total_minutes" column)
MINUTE = 4320

# Attribute to plot (the column prefix without the "_<r_in>_to_<r_out>" suffix).
# Examples: num_cells_radius, num_tumor_cells_radius, tumor_cells_type1_radius,
# ..., tumor_cells_type4_radius, tumor_cells_type5_dead_radius,
# dead_tumor_cells_lack_of_oxygen_radius, dead_tumor_cells_lack_of_glucose_radius,
# dead_tumor_cells_random_natural_causes_radius,
# dead_tumor_cells_cartcell_kill_radius, num_alive_cart_radius,
# num_dead_cart_radius, average_oncoprotein_radius,
# average_oxygen_cancer_cells_radius, average_oxygen_all_cells_radius,
# average_glucose_cancer_cells_radius, average_glucose_all_cells_radius
ATTRIBUTE = "tumor_cells_type5_dead_radius"

# Divide the values by this factor to change units
# (e.g. 585 converts oxygen from mmHg to mol/m3)
DIVIDING_FACTOR = 1

# True: divide counts by the ring area, so the value is expressed per um^2
# instead of as the total amount inside the ring
APPLY_RADIAL_NORMALIZATION = True

# Y axis limits (None = automatic)
Y_MIN = 0
Y_MAX = None


def load_profile(csv_path, attribute, minute):
    """Returns (radii, inner radii, outer radii, values) sorted by radius."""
    df = pd.read_csv(csv_path)
    rows = df[df["total_minutes"] == minute]
    if rows.empty:
        raise SystemExit(f"Minute {minute} not found in {csv_path}")
    row = rows.iloc[0]

    pattern = re.compile(rf"^{re.escape(attribute)}_(\d+)_to_(\d+)$")
    points = []
    for col in df.columns:
        match = pattern.match(col)
        if match and pd.notna(row[col]):
            r_in, r_out = int(match.group(1)), int(match.group(2))
            points.append(((r_in + r_out) / 2, r_in, r_out, row[col]))

    if not points:
        raise SystemExit(f"No columns found matching attribute '{attribute}'")

    points.sort(key=lambda p: p[0])
    radii, r_in, r_out, values = (np.array(c, dtype=float) for c in zip(*points))
    return radii, r_in, r_out, values


def main():
    radii, r_in, r_out, values = load_profile(CSV_PATH, ATTRIBUTE, MINUTE)

    if APPLY_RADIAL_NORMALIZATION:
        values = values / (np.pi * (r_out**2 - r_in**2))
    values = values / DIVIDING_FACTOR

    print(f"Attribute: {ATTRIBUTE}")
    print(f"Minute: {MINUTE}")
    print(f"Radius range with data: {radii[0]:g} to {radii[-1]:g}")
    print(f"Number of points: {len(radii)}")
    print(f"Radial normalization applied: {APPLY_RADIAL_NORMALIZATION}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(radii, values, marker="o")
    ax.set_xlabel("Radius (um)")
    ax.set_ylabel(ATTRIBUTE + (" (per um^2)" if APPLY_RADIAL_NORMALIZATION else ""))
    ax.set_title(f"{ATTRIBUTE} vs radius - minute {MINUTE}")
    ax.grid(True, alpha=0.3)
    if Y_MIN is not None or Y_MAX is not None:
        ax.set_ylim(bottom=Y_MIN, top=Y_MAX)
    fig.tight_layout()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / f"{ATTRIBUTE}_vs_radius_minute_{MINUTE}.png"
    fig.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
