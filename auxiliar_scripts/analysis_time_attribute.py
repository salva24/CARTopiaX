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

This file plots one attribute of the simulation as a function of time, using
output/final_data.csv (one row per saved time step, one column per attribute).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
CSV_PATH = BASE_DIR / "output" / "final_data.csv"
RESULTS_DIR = Path(__file__).resolve().parent / "out"  # destination folder

# Attribute to plot (a column of final_data.csv). Options: tumor_radius,
# num_cells, num_tumor_cells, tumor_cells_type1, ..., tumor_cells_type4,
# tumor_cells_type5_dead, dead_tumor_cells_lack_of_oxygen,
# dead_tumor_cells_lack_of_glucose, dead_tumor_cells_random_natural_causes,
# dead_tumor_cells_cartcell_kill, num_alive_cart, num_dead_cart,
# average_oncoprotein, average_oxygen_cancer_cells, average_oxygen_all_cells,
# average_glucose_cancer_cells, average_glucose_all_cells,
# average_radius_distance_living_cart_cells
ATTRIBUTE = "average_oncoprotein"

# Time unit of the x axis: "days", "hours" or "minutes"
TIME_UNIT = "days"

# Divide the values by this factor to change units
# (e.g. 585 converts oxygen from mmHg to mol/m3)
DIVIDING_FACTOR = 1

# Y axis limits (None = automatic)
Y_MIN = None
Y_MAX = None

COLOR = "#e41a1c"


def main():
    df = pd.read_csv(CSV_PATH)

    time_col = f"total_{TIME_UNIT}"
    if time_col not in df.columns:
        raise SystemExit(f"Invalid TIME_UNIT '{TIME_UNIT}' (use days, hours or minutes)")
    if ATTRIBUTE not in df.columns:
        raise SystemExit(
            f"Attribute '{ATTRIBUTE}' not found in {CSV_PATH}. "
            f"Available: {', '.join(df.columns[3:])}"
        )

    time_points = df[time_col]
    values = df[ATTRIBUTE] / DIVIDING_FACTOR

    print(f"Attribute: {ATTRIBUTE}")
    print(f"Time range: {time_points.iloc[0]:g} to {time_points.iloc[-1]:g} {TIME_UNIT}")
    print(f"Number of points: {len(df)}")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time_points, values, color=COLOR, linestyle="-", linewidth=2)
    ax.set_xlabel(f"Time ({TIME_UNIT})")
    ax.set_ylabel(ATTRIBUTE)
    ax.set_title(f"{ATTRIBUTE} over time")
    ax.grid(True, linestyle="--", alpha=0.7)
    if Y_MIN is not None or Y_MAX is not None:
        ax.set_ylim(bottom=Y_MIN, top=Y_MAX)
    fig.tight_layout()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / f"{ATTRIBUTE}_vs_time.png"
    fig.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
