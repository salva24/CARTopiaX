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

This file contains the Bayesian optimization workflow used to calibrate
the model parameters, developed under Google Summer of Code (GSoC)
for the compiler-research.org organization.
"""

import json
import shutil
import subprocess
from pathlib import Path

# --- Adjust these paths/values to your environment ---
# Directory for the visualization configuration file bdm.toml
BDM_TOML_PATH = Path("./bdm.toml")
BASE_DIR = Path(__file__).resolve().parent.parent
# Directory for the parameters file params.json
PARAMS_PATH = BASE_DIR / "params.json"

BIODYNAMO_DIR = " /home/usuario/Desktop/biodynamo/build/bin/thisbdm.sh"
NUMBER_EXECUTIONS = 10

SIM_DIR = BASE_DIR / "output"        # folder where the ABM writes its results
DRAFT_DIR = BASE_DIR / "draft"       # base destination folder


def run_ABM(seed):

    config={
  "seed": seed,
  "bound_space_toplogy": "closed",
  "num_radius_intervals": 20,
  "lateral_oxygen_production_min_z": -300.0,
  "lateral_oxygen_production_max_z": 300.0,
  "min_initial_z_substances": -300.0,
  "max_initial_z_substances": 300.0,
  "diffuse_oxygen_on_z_axis": False,
  "diffuse_glucose_on_z_axis": False,
  "output_performance_statistics": False,
  "total_minutes_to_simulate": 4320,
  "output_csv_interval": 600,
  "bounded_space_length": 6500.0,
  "tumor_shape": "cylinder",
  "output_information_dependent_on_radius": True,
  "cylindrical_tumor_radius": 3000.0,
  "max_radius_analysis_csv_dependent_on_radius": 3000.0,
  "cylindrical_tumor_height": 100.0,
  "initial_number_of_cylindrical_tumor_cells": 28000,
  "default_volume_new_tumor_cell": 1468.0,
  "std_volume_new_tumor_cell": 166,
  "min_volume_new_tumor_cell": 1136.0,
  "max_volume_new_tumor_cell": 1800.0,
  "default_volume_new_cart_cell": 269.0,
  "std_volume_new_cart_cell": 33.625,
  "min_volume_new_cart_cell": 201.75,
  "max_volume_new_cart_cell": 336.25,
  "oncoprotein_mean": 1.0,
  "oncoprotein_standard_deviation": 0.0,
  "initial_oxygen_level": 0.0,
  "oxygen_reference_level": 165.0,
  "default_oxygen_consumption_tumor_cell": 77.47,
  "default_glucose_consumption_tumor_cell": 0.045,
  "diffusion_coefficient_oxygen": 180000.0,
  "decay_constant_oxygen": 0.01,
  "time_apoptosis": 6000.0,
  "time_lysis": 6000.0,
  "treatment": {
    "0": 0
  },
  "average_time_transformation_random_rate": 72,
  "standard_deviation_transformation_random_rate": 15.0,
  "decay_constant_glucose": 0.0005,
  "oxygen_saturation_for_proliferation": 13.74,
  "oxygen_limit_for_proliferation": 5.9,
  "oxygen_limit_for_necrosis_maximum": 0.0,
  "oxygen_limit_for_necrosis": 45.0,
  "maximum_necrosis_lack_of_oxygen_rate": 0.0000216,
  "reduction_consumption_dead_cells": 0.0,
  "basal_death_probability_cancer_cells": 0.000005,
  "bounded_space_min_allowed_z": -50.0,
  "bounded_space_max_allowed_z": 50.0,
  "bounded_space_max_allowed_radius": 3250.0,
  "minimum_distance_from_tumor_to_spawn_cart": 0.0,
  "add_immunostimulatory_factor": True,
  "add_glucose": True,
  "diffusion_coefficient_glucose": 7800,
  "initial_glucose_level": 24.98,
  "max_radius_glucose_initialization": 3100.0,  
  "default_oxygen_consumption_cart": 7,
  "default_glucose_consumption_cart": 0.07,
  "glucose_saturation_for_tumor_cell_growth": 24.98,
  "glucose_limit_for_tumor_cell_growth": 0,
  "glucose_limit_for_death": 5,
  "glucose_limit_for_death_maximum": 0,
  "maximum_death_lack_of_glucose_rate": 0.000175,
  "minimum_tumor_cell_target_volume_fraction_for_division": 0.9,
  "avg_migration_bias_cart": 0,
  "std_migration_bias_cart": 0.135,
  "persistence_time_cart": 0,
  "kill_rate_cart": 0.521,
  "adhesion_rate_cart": 0.00032,
  "diffusion_coefficient_immunostimulatory_factor": 100000,
  "decay_constant_immunostimulatory_factor": 0.00001
}

    # Save the config parameters for this run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)

    # Load the BioDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(
        ["bash", "-c", f"source {BIODYNAMO_DIR} && bdm run"], check=True
    )


def copy_output_to_draft(seed):
    """Copy the contents of SIM_DIR to draft/execution_seed_{seed}."""
    dest_dir = DRAFT_DIR / f"execution_seed_{seed}"

    if not SIM_DIR.exists():
        raise FileNotFoundError(f"Results folder not found: {SIM_DIR}")

    if dest_dir.exists():
        shutil.rmtree(dest_dir)  # avoid mixing results from previous runs

    shutil.copytree(SIM_DIR, dest_dir)
    print(f"[seed={seed}] Results copied to: {dest_dir}")


def run_multiple_seeds(N):
    """Run the ABM with seeds 0..N-1 and save each result separately."""
    DRAFT_DIR.mkdir(parents=True, exist_ok=True)

    # for seed in range(N):
    for seed in range(N):
        print(f"--- Running ABM with seed={seed} ---")
        run_ABM(seed)
        copy_output_to_draft(seed)

    print(f"Done: {N} runs completed (seeds 0 to {N - 1}).")

# Function to set the export value in the bdm.toml file
def set_export_value(path: Path, value: bool):
    lines = path.read_text().splitlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("export"):
            new_lines.append(f"export = {str(value).lower()}")
        else:
            new_lines.append(line)

    path.write_text("\n".join(new_lines))

if __name__ == "__main__":

    # Save the original content of the bdm.toml and params.json files to restore them later
    original_bdm = BDM_TOML_PATH.read_text()
    original_params = PARAMS_PATH.read_text()

    try:
        # Change the export value in the bdm.toml file to False to disable visualization
        set_export_value(BDM_TOML_PATH, False)

        run_multiple_seeds(N=NUMBER_EXECUTIONS)

    finally:
        # Restore the original content of the bdm.toml and params.json files always, even if an error occurs during the optimization process
        BDM_TOML_PATH.write_text(original_bdm)
        PARAMS_PATH.write_text(original_params)