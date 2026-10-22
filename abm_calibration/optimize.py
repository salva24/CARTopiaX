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
import logging
import subprocess
from pathlib import Path

import numpy as np
import optuna
import pandas as pd

# Change this: File Parameters for the desired experiment
EXPERIMENT_ID = 0
# "total": compares against target_data/final_data.csv
# "radius": compares against target_data/data_dependent_on_radius_tumor.csv
# "custom": implement your own error function in compute_custom_error()
MODE = "total"

SEED = 42
# You can set the number of trials to 0 to skip the optimization and just load the best result from the database
NUMBER_OF_TRIALS = 8
# Number of Monte Carlo simulations to run for each trial. Use one for an aproximation of the error with a single montecarlo run
NUMBER_MONTE_CARLO = 3
# BioDynaMo directory to execute the comand source thisbdm.sh, you can change it to your own path
BIODYNAMO_DIR = "/home/usuario/Desktop/biodynamo/build/bin/thisbdm.sh"

# Other Hyperparameters for the optimization process
# LOGGING
logging.basicConfig(level=logging.INFO)

# Directory for the visualization configuration file bdm.toml
BDM_TOML_PATH = Path("./bdm.toml")
# Directory for the parameters file params.json
PARAMS_PATH = Path(__file__).resolve().parent.parent / "params.json"

# Directories with the simulation output and the target data
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"
TARGET_DIR = Path(__file__).resolve().parent / "target_data"

# Directory for this experiment
EXPERIMENT_DIR = Path("abm_calibration") / f"experiment_{EXPERIMENT_ID}"
EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)

#######################################################
# Change this: Simulations to run the ABM with the given parameters
#######################################################


# Function to run the ABM with the given parameters
def run_ABM(params, seed):
    # Change this: parameter to be optimized in the ABM simulation
    initial_oxygen_level = params["initial_oxygen_level"]
    default_oxygen_consumption_tumor_cell = params[
        "default_oxygen_consumption_tumor_cell"
    ]

    # Change this configuration for the ABM run
    # (for MODE = "radius" it must set output_information_dependent_on_radius: True)
    config = {
        "seed": seed,
        "output_performance_statistics": False,
        "total_minutes_to_simulate": 1440,
        "initial_tumor_radius": 40.0,
        "treatment": {"0": 50},
        "initial_oxygen_level": initial_oxygen_level,
        "default_oxygen_consumption_tumor_cell": default_oxygen_consumption_tumor_cell,
    }

    # Save the config parameters for the run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)

    # Load the ByoDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(["bash", "-c", f"source {BIODYNAMO_DIR} && bdm run"], check=True)


#######################################################
# Change this: Error functions (examples)
#######################################################


# Error using the data of the whole tumor as a function of time (final_data.csv)
def compute_error_total():
    metric = "average_oxygen_all_cells"  # Change this: choose the metric to compare, e.g., "average_oxygen_all_cells", "tumor_radius", etc.
    return mse_against_target("final_data.csv", [metric])


# Error using the data as a function of the radius at a fixed time (data_dependent_on_radius_tumor.csv)
def compute_error_radius():
    csv_name = "data_dependent_on_radius_tumor.csv"
    metric = "tumor_cells_type5_dead"  # Change this: prefix of the metric columns, e.g., "tumor_cells_type5_dead"
    minute = 4320  # Change this: fixed time (total_minutes) at which the radial profile is compared

    target = TARGET_DIR / csv_name
    if not target.exists():
        logging.error("Missing CSV: %s", target)
        return float("inf")

    header = pd.read_csv(target, nrows=0).columns
    columns = [c for c in header if c.startswith(f"{metric}_radius_")]
    if not columns:
        logging.error("No columns for metric %s in %s", metric, target)
        return float("inf")
    return mse_against_target(csv_name, columns, minute=minute)


# Change this: write your own error. Example: dead cells at the center of a tumor
# of 3000 micrometers of radius split into 20 measurement rings at minute 4320 (72 hours)
def compute_custom_error():
    sim = OUTPUT_DIR / "data_dependent_on_radius_tumor.csv"
    if not sim.exists():
        logging.error("Missing simulation CSV: %s", sim)
        return float("inf")

    df_s = pd.read_csv(
        sim, usecols=["total_minutes", "tumor_cells_type5_dead_radius_0_to_150"]
    )
    row = df_s[df_s["total_minutes"] == 4320].iloc[0]
    target_value_center = 21
    return abs(
        float(row["tumor_cells_type5_dead_radius_0_to_150"]) - target_value_center
    )


#######################################################
# Change this: Objective function for the Optuna optimization process
#######################################################


def objective(trial):
    # Change this: Define the parameters to be optimized and their ranges
    params = {
        "initial_oxygen_level": trial.suggest_float("initial_oxygen_level", 30, 40),
        "default_oxygen_consumption_tumor_cell": trial.suggest_float(
            "default_oxygen_consumption_tumor_cell", 7, 14
        ),
    }

    logging.info(f"Trial {trial.number} | params={params}")

    # Compute the error as the average of the errors from multiple Monte Carlo simulations varying the seed
    total_error = 0
    for seed in np.random.randint(0, 10000, NUMBER_MONTE_CARLO):
        run_ABM(params, int(seed))

        error = compute_error()
        total_error += error

    error = total_error / NUMBER_MONTE_CARLO

    logging.info(f"Trial {trial.number} | error={error}")

    return error


#######################################################
# Auxiliary functions
#######################################################


# Function to compute the error between the ABM simulation results and the target data
def compute_error():
    if MODE == "total":
        return compute_error_total()
    if MODE == "radius":
        return compute_error_radius()
    if MODE == "custom":
        return compute_custom_error()
    raise ValueError(f"Unknown MODE: {MODE}")


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


# MSE between simulation and target for the given csv file and columns, merged on total_minutes
# (if minute is given, only that time is compared)
def mse_against_target(csv_name, columns, minute=None):
    target = TARGET_DIR / csv_name
    sim = OUTPUT_DIR / csv_name

    for path in (target, sim):
        if not path.exists():
            logging.error("Missing CSV: %s", path)
            return float("inf")

    df_t = pd.read_csv(target, usecols=["total_minutes", *columns])
    df_s = pd.read_csv(sim, usecols=["total_minutes", *columns])

    if minute is not None:
        df_t = df_t[df_t["total_minutes"] == minute]
        df_s = df_s[df_s["total_minutes"] == minute]

    merged = pd.merge(
        df_t, df_s, on="total_minutes", how="inner", suffixes=("_t", "_s")
    )
    if merged.empty:
        logging.error("No common minutes between target and simulation")
        return float("inf")

    squared = [(merged[f"{c}_s"] - merged[f"{c}_t"]) ** 2 for c in columns]
    return float(pd.concat(squared, axis=1).stack().mean())


#######################################################
# Main
#######################################################

if __name__ == "__main__":
    # Fix the random seed for reproducibility
    np.random.seed(SEED)

    # Create an Optuna study to optimize the parameters of the ABM
    study = optuna.create_study(
        study_name="abm_calibration",
        storage=f"sqlite:///{EXPERIMENT_DIR / 'abm_optuna.db'}",
        load_if_exists=True,
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
    )

    # Save the original content of the bdm.toml and params.json files to restore them later
    original_bdm = BDM_TOML_PATH.read_text()
    original_params = PARAMS_PATH.read_text()

    try:
        # Change the export value in the bdm.toml file to False to disable visualization
        set_export_value(BDM_TOML_PATH, False)

        study.optimize(objective, n_trials=NUMBER_OF_TRIALS)

        print("\nBEST RESULT")
        print("Value:", study.best_value)
        print("Params:", study.best_params)

        df = study.trials_dataframe()
        df.to_csv(EXPERIMENT_DIR / "optuna_results.csv", index=False, mode="w")

    finally:
        # Restore the original content of the bdm.toml and params.json files always, even if an error occurs during the optimization process
        BDM_TOML_PATH.write_text(original_bdm)
        PARAMS_PATH.write_text(original_params)
