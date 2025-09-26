/*
 * Copyright 2025 compiler-research.org, Salvador de la Torre Gonzalez, Luciana
 * Melina Luque
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *     SPDX-License-Identifier: Apache-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * This file contains a model developed under Google Summer of Code (GSoC)
 * for the compiler-research.org organization.
 */

#ifndef TUMOR_HYPERPARAMS_H_
#define TUMOR_HYPERPARAMS_H_

#include "core/param/param_group.h"
#include "core/real_t.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>
#include <map>
#include <nlohmann/json.hpp>

namespace bdm {

///
/// Some constants defined for code legibility
///
/// Twice Pi
constexpr real_t kTwicePi = 2. * Math::kPi;
/// Constant for division by 2.0
constexpr real_t kHalf = 2.0;
/// Epsilon for avoiding division by 0
constexpr real_t kEpsilon = 1e-10;
/// kEpsilon distance
constexpr real_t kEpsilonDistance = 1e-5;
/// Large time to avoid division by 0
constexpr real_t kTimeTooLarge = 1e100;
/// Minutes in an hour
constexpr real_t kMinutesInAnHour = 60.0;
/// Hours in a day
constexpr real_t kHoursInADay = 24.0;

/// Contains the default values of the hyperparameters used in the simulation.
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

    ///
    /// General Hyperparameters
    ///

    /// Seed for random number generation
    int kSeed = 1;
    /// Output Performance Statistics
    bool kOutputPerformanceStatistics = false;
    /// Total simulation time in minutes (30 days)
    int kTotalMinutesToSimulate = 43200;
    /// Initial radius of the spherical tumor (group of cancer cells) in micrometers
    real_t kInitialRadiusTumor = 150;
    /// Length of the bounded space in micrometers
    int kBoundedSpaceLength = 1000;

    /// Treatment Dosages
    ///
    /// Specifies the CAR-T cell infusion schedule as a map where:
    ///   - The key represents the day of treatment (starting from day 0).
    ///   - The value represents the number of CAR-T cells administered on that day.
    /// Example: On day 0 and day 8, 3957 CAR-T cells are introduced (matching the
    /// initial tumor cell count). {0, 3957}, {8, 3957}}
    // avoid cppcoreguidelines-avoid-magic-numbers warning
    std::map<int, int> kTreatment = {{0, 3957}, {8, 3957}};

    /// Time steps
    /// Time step for substances secretion/consumption in minutes
    real_t kDtSubstances = 0.01;
    /// Time step for the cell mechanics in minutes
    real_t kDtMechanics = 0.1;
    /// Time step for the cell cycle in minutes
    real_t kDtCycle = 6;
    /// General time step for the simulation: it is the same as kDtMechanics
    real_t kDt = 0.1;

    /// Output little summary each half kOutputCsvInterval.
    /// By default half a day (12h * 60min / kDt)
    int kOutputCsvInterval = 7200;

    /// Apoptotic cells volume change
    real_t kVolumeRelaxationRateCytoplasmApoptotic = 0.0166667;
    real_t kVolumeRelaxationRateNucleusApoptotic = 0.00583333;
    real_t kVolumeRelaxationRateFluidApoptotic = 0.0;
    /// Time in minutes until an apoptotic cell is removed from the simulation
    real_t kTimeApoptosis = 516;
    /// Reduction of consumption rate of dead cells when they enter necrosis
    real_t kReductionConsumptionDeadCells = 0.1;

    /// Chemicals
    ///  Number of voxels per axis for the substances grid
    int kResolutionGridSubstances = 50;
    /// Diffusion coefficient of oxygen in μm²/min
    real_t kDiffusionCoefficientOxygen = 100000;
    /// Decay constant of oxygen in min⁻¹
    real_t kDecayConstantOxygen = 0.1;
    /// Diffusion coefficient of immunostimulatory factor in μm²/min
    real_t kDiffusionCoefficientImmunostimulatoryFactor = 1000;
    /// Decay constant of immunostimulatory factor in min⁻¹
    real_t kDecayConstantImmunostimulatoryFactor = 0.016;
    /// Reference level of oxygen at the boundaries in mmHg
    real_t kOxygenReferenceLevel = 38;
    /// Initial oxygen concentration in each voxel in mmHg
    real_t kInitialOxygenLevel = 38;
    /// Oxygen saturation in the microenvironment in mmHg
    real_t kOxygenSaturation = 30;

    /// Forces
    ///  Repulsion coeficient between tumor cells
    real_t kRepulsionTumorTumor = 10;
    /// Repulsion coeficient between CAR-T cells
    real_t kRepulsionCartCart = 50;
    /// Repulsion coeficient between CAR-T cells and tumor cells
    real_t kRepulsionCartTumor = 50;
    /// Repulsion coeficient between tumor cells and CAR-T cells
    real_t kRepulsionTumorCart = 10;
    /// Maximum relative adhesion distance for cell interactions
    real_t kMaxRelativeAdhesionDistance = 1.25;
    /// Adhesion coeficient between tumor cells
    real_t kAdhesionTumorTumor = 0.4;
    /// Adhesion coeficient between CAR-T cells
    real_t kAdhesionCartCart = 0;
    /// Adhesion coeficient between CAR-T cells and tumor cells
    real_t kAdhesionCartTumor = 0;
    /// Adhesion coeficient between tumor cells and CAR-T cells
    real_t kAdhesionTumorCart = 0;
    /// Box length for mechanics calculations (in micrometers). Smaller values
    /// improve simulation efficiency
    real_t kLengthBoxMechanics = 22;
    // Coefficientes for the two step Adams-Bashforth approximation of the time
    // derivative for position position(t + dt) ≈ position(t) + dt * [ 1.5 *
    // velocity(t) - 0.5 * velocity(t - dt) ]
    /// Coefficient for the current time step in the Adams-Bashforth method (dt
    /// * 1.5)
    real_t kDnew = 0.15;
    /// Coefficient for the previous time step in the Adams-Bashforth method (dt *
    /// -0.5)
    real_t kDold = -0.05;
    /// maximum squared distance to avoid CAR-T pushing
    /// tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
    /// CAR-T cell will only move to the tumor cell with the adhesion forces
    /// (radiusCART + radiusTumorCell + 0.1 to avoid numerical errors)**2
    real_t kMaxSquaredDistanceCartMovingTowardsTumorCell = 317.746;

    ///
    /// TumorCell Hyperparameters
    ///

    /// Rate of secretion of immunostimulatory factor of tumor cells per minute
    real_t kRateSecretionImmunostimulatoryFactor = 10;
    /// Saturation density of immunostimulatory factor for tumor cells
    real_t kSaturationDensityImmunostimulatoryFactor = 1;
    /// Mean level of oncoprotein expression in tumor cells
    real_t kOncoproteinMean = 1;
    /// Standard deviation of oncoprotein expression in tumor cells
    real_t kOncoproteinStandardDeviation = 0.25;
    /// Oxygen saturation level in tumor cells for proliferation
    real_t kOxygenSaturationInProliferation = 38;
    /// Limit of oxygen level for tumor cell proliferation
    real_t kOxygenLimitForProliferation = 10;
    /// Limit of oxygen to start causing necrosis
    real_t kOxygenLimitForNecrosis = 5;
    /// Limit of oxygen to maximum necrosis probability
    real_t kOxygenLimitForNecrosisMaximum = 2.5;
    /// Time in minutes until a lysed necrotic cell is removed from the simulation
    real_t kTimeLysis = 86400;
    /// Maximum rate per minute of necrosis for tumor cells in case of hypoxia with
    /// 0 oxygen
    real_t kMaximumNecrosisRate = 0.00277778;
    /// Default oxygen consumption rate of tumor cell
    real_t kDefaultOxygenConsumptionTumorCell = 10;
    /// Volume parameters
    /// Default total volume of a new tumor cell in μm³
    real_t kDefaultVolumeNewTumorCell = 2494;
    /// Default volume of the nucleus of a new tumor cell in μm³
    real_t kDefaultVolumeNucleusTumorCell = 540;
    /// Default fraction of fluid volume in a new tumor cell
    real_t kDefaultFractionFluidTumorCell = 0.75;
    /// Average time for transformation Random Rate in hours
    real_t kAverageTimeTransformationRandomRate = 38.6;
    /// Standard Deviation for transformation Random Rate in hours
    real_t kStandardDeviationTransformationRandomRate = 3.7;
    /// Average adhesion time in minutes for Tumor Cell under CAR-T attack before
    /// escaping
    real_t kAdhesionTime = 60;
    /// Min oncoprotein level to be killed by a CAR-T cell
    real_t kOncoproteinLimit = 0.5;
    /// Max oncoprotein level
    real_t kOncoproteinSaturation = 2.0;
    /// Do not modify this line: difference between saturation and limit
    real_t kOncoproteinDifference = 1.5;

    /// volume relaxation rate (min^-1) for each state
    real_t kVolumeRelaxationRateAliveCytoplasm = 0.00216667;
    real_t kVolumeRelaxationRateAliveNucleus = 0.00366667;
    real_t kVolumeRelaxationRateAliveFluid = 0.0216667;

    real_t kVolumeRelaxationRateCytoplasmNecroticSwelling = 5.33333e-05;
    real_t kVolumeRelaxationRateNucleusNecroticSwelling = 0.000216667;
    real_t kVolumeRelaxationRateFluidNecroticSwelling = 0.000833333;

    real_t kVolumeRelaxationRateCytoplasmNecroticLysed = 5.33333e-05;
    real_t kVolumeRelaxationRateNucleusNecroticLysed = 0.000216667;
    real_t kVolumeRelaxationRateFluidNecroticLysed = 0.000833333;

    /// Thresholds in oncoprotein levels for differentiating 4 cancer cell types
    real_t kThresholdCancerCellType1 = 1.5;
    real_t kThresholdCancerCellType2 = 1.0;
    real_t kThresholdCancerCellType3 = 0.5;
    real_t kThresholdCancerCellType4 = 0.0;

    ///
    /// CAR-T Cell Hyperparameters
    ///

    /// Average time in minutes until a CAR-T cell dies
    real_t kAverageMaximumTimeUntillApoptosisCart = 86400;
    /// Default oxygen consumption rate of CAR-T cell
    real_t kDefaultOxygenConsumptionCarT = 1;
    /// Volume parameters
    ///  Default total volume of a new CAR-T cell in μm³
    real_t kDefaultVolumeNewCarTCell = 2494;

    /// How often a CAR-T cell tries to kill an attached cancer cell in 1/min
    real_t kKillRateCart = 0.06667; 
    /// How often a CAR-T cell tries to attach to a cancer cell in 1/min
    real_t kAdhesionRateCart = 0.013;
    /// Maximum adhesion distance between CAR-T and tumor cells in micrometers
    real_t kMaxAdhesionDistanceCart = 18;
    /// Minimum adhesion distance between CAR-T and tumor cells in micrometers
    real_t kMinAdhesionDistanceCart = 14;

    /// Minimum distance in micrometers from the tumor for spawning CAR-T cells
    real_t kMinimumDistanceCarTFromTumor = 50;

    /// Motility parameters
    /// Average persistence time before CAR-T cell moves
    real_t kPersistenceTimeCart = 10;  // 10 minutes
    /// Higher bias (\in [0,1]) makes CAR-T movement more directed toward
    /// immunostimulatory factor source; while a bias of 0 makes the movement random
    real_t kMigrationBiasCart = 0.5;
    /// Migration speed
    real_t kMigrationSpeedCart = 5;
    /// Elastic constant
    real_t kElasticConstantCart = 0.01;

    ///
    /// Fixed Constants that should not be directly changed
    ///
    
    /// Number of steps per cycle step Needs to be computed to avoid errors with
    /// fmod
    int kStepsPerCycle = 60;
    /// Steps in one day
    size_t kStepsOneDay = 14400;
    /// Volume of a single voxel in μm³
    real_t kVoxelVolume = 8000;
    /// 1-kMigrationBiasCart
    real_t kMigrationOneMinusBiasCart = 0.5;
    /// Probability of a CAR-T cell to migrate in a given
    /// mechanical time step
    real_t kMotilityProbability = 0.01;
    /// Probability of a Tumor cell to escape in a given
    /// mechanical time step
    real_t kProbabilityEscape = 0.00166667;

    /// Maximum adhesion distance squared
    real_t kSquaredMaxAdhesionDistanceCart = 324;
    /// Difference between min and max adhesion distance
    real_t kDifferenceCartAdhesionDistances = 4;

    /// Radius tumor cell
    real_t kRadiusTumorCell = 8.41271;

    /// Radius cart cell
    real_t kRadiusCarTCell = 8.41271;

    /// Max Distance for considering two cells as neighbours for force calculations
    /// in μm (twice cell radius times kMaxRelativeAdhesionDistance + 0.1 to avoid
    /// mismatch because of numerical errors)**2
    real_t kSquaredMaxDistanceNeighborsForce = 446.552;

    
    ///
    /// Function to load the parameters
    ///

    /// Function to read a JSON file and load the hyperparameters in it.
    /// The parameters that are not found in the file will be assigned their default value
    /// It is important that this function is called before starting the simulation
    void LoadParams(const std::string& filename);

    ///
    /// Function to print all the parameters with their values
    ///

    /// Print all the parameters with their values
    void PrintParams() const;
};




}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_
