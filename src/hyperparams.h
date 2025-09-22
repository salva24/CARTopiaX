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

#include "core/real_t.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>
#include <map>

namespace bdm {

/// This file contains the default values of the hyperparameters used in the
/// simulation.

///
/// General Hyperparameters
///

/// Seed for random number generation
constexpr int kSeed = 1;
/// Output Performance Statistics
constexpr bool kOutputPerformanceStatistics = true;
/// Total simulation time in minutes (30 days)
constexpr int kTotalMinutesToSimulate = 30 * 24 * 60;
/// Initial radius of the spherical tumor (group of cancer cells) in micrometers
constexpr real_t kInitialRadiusTumor = 150;
/// Length of the bounded space in micrometers
constexpr int kBoundedSpaceLength = 1000;

/// Treatment Dosages
///
/// Specifies the CAR-T cell infusion schedule as a map where:
///   - The key represents the day of treatment (starting from day 0).
///   - The value represents the number of CAR-T cells administered on that day.
/// Example: On day 0 and day 8, 3957 CAR-T cells are introduced (matching the
/// initial tumor cell count).
// avoid cppcoreguidelines-avoid-magic-numbers warning
inline const std::map<size_t, size_t> gKTreatment = {
    {0, 3957},
    {8,
     3957}};  // NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

/// Time steps
/// Time step for substances secretion/consumption in minutes
constexpr real_t kDtSubstances = 0.01;
/// Time step for the cell mechanics in minutes
constexpr real_t kDtMechanics = 0.1;
/// Time step for the cell cycle in minutes
constexpr real_t kDtCycle = 6.0;
/// General time step for the simulation: it is the same as kDtMechanics
constexpr real_t kDt = kDtMechanics;

/// Output little summary each half a day
constexpr int kOutputCsvInterval = 12 * 60 / kDt;

/// Apoptotic cells volume change
constexpr real_t kVolumeRelaxationRateCytoplasmApoptotic = 1.0 / 60.0;
constexpr real_t kVolumeRelaxationRateNucleusApoptotic = 0.35 / 60.0;
constexpr real_t kVolumeRelaxationRateFluidApoptotic = 0.0;
/// Time in minutes until an apoptotic cell is removed from the simulation
constexpr real_t kTimeApoptosis = 8.6 * 60;
/// Reduction of consumption rate of dead cells when they enter necrosis
constexpr real_t kReductionConsumptionDeadCells = 0.1;

/// Chemicals
///  Number of voxels per axis for the substances grid
constexpr int kResolutionGridSubstances = 50;
/// Diffusion coefficient of oxygen in μm²/min
constexpr real_t kDiffusionCoefficientOxygen = 100000;
/// Decay constant of oxygen in min⁻¹
constexpr real_t kDecayConstantOxygen = 0.1;
/// Diffusion coefficient of immunostimulatory factor in μm²/min
constexpr real_t kDiffusionCoefficientImmunostimulatoryFactor = 1000;
/// Decay constant of immunostimulatory factor in min⁻¹
constexpr real_t kDecayConstantImmunostimulatoryFactor = 0.016;
/// Reference level of oxygen at the boundaries in mmHg
constexpr real_t kOxygenReferenceLevel = 38.;
/// Initial oxygen concentration in each voxel in mmHg
constexpr real_t kInitialOxygenLevel = 38.0;
/// Oxygen saturation in the microenvironment in mmHg
constexpr real_t kOxygenSaturation = 30.0;

/// Forces
///  Repulsion coeficient between tumor cells
constexpr real_t kRepulsionTumorTumor = 10.0;
/// Repulsion coeficient between CAR-T cells
constexpr real_t kRepulsionCartCart = 50.0;
/// Repulsion coeficient between CAR-T cells and tumor cells
constexpr real_t kRepulsionCartTumor = 50.0;
/// Repulsion coeficient between tumor cells and CAR-T cells
constexpr real_t kRepulsionTumorCart = 10.0;
/// Maximum relative adhesion distance for cell interactions
constexpr real_t kMaxRelativeAdhesionDistance = 1.25;
/// Adhesion coeficient between tumor cells
constexpr real_t kAdhesionTumorTumor = 0.4;
/// Adhesion coeficient between CAR-T cells
constexpr real_t kAdhesionCartCart = 0.0;
/// Adhesion coeficient between CAR-T cells and tumor cells
constexpr real_t kAdhesionCartTumor = 0.0;
/// Adhesion coeficient between tumor cells and CAR-T cells
constexpr real_t kAdhesionTumorCart = 0.0;
/// Box length for mechanics calculations (in micrometers). Smaller values
/// improve simulation efficiency
const real_t gKLengthBoxMechanics = 22;
// Coefficientes for the two step Adams-Bashforth approximation of the time
// derivative for position position(t + dt) ≈ position(t) + dt * [ 1.5 *
// velocity(t) - 0.5 * velocity(t - dt) ]
/// Coefficient for the current time step in the Adams-Bashforth method (dt
/// * 1.5)
constexpr real_t kDnew = 1.5 * kDtMechanics;
/// Coefficient for the previous time step in the Adams-Bashforth method (dt *
/// -0.5)
constexpr real_t kDold = -0.5 * kDtMechanics;

///
/// TumorCell Hyperparameters
///

/// Rate of secretion of immunostimulatory factor of tumor cells per minute
constexpr real_t kRateSecretionImmunostimulatoryFactor = 10.0;
/// Saturation density of immunostimulatory factor for tumor cells
constexpr real_t kSaturationDensityImmunostimulatoryFactor = 1.0;
/// Mean level of oncoprotein expression in tumor cells
constexpr real_t kOncoproteinMean = 1.0;
/// Standard deviation of oncoprotein expression in tumor cells
constexpr real_t kOncoproteinStandardDeviation = 0.25;
/// Oxygen saturation level in tumor cells for proliferation
constexpr real_t kOxygenSaturationInProliferation = 38.0;
/// Limit of oxygen level for tumor cell proliferation
constexpr real_t kOxygenLimitForProliferation = 10.0;
/// Limit of oxygen to start causing necrosis
constexpr real_t kOxygenLimitForNecrosis = 5.0;
/// Limit of oxygen to maximum necrosis probability
constexpr real_t kOxygenLimitForNecrosisMaximum = 2.5;
/// Time in minutes until a lysed necrotic cell is removed from the simulation
constexpr real_t kTimeLysis = 60 * 24 * 60.;
/// Rate of cell division in 1/min
constexpr real_t kDivisionRate = 0.02717 / 60.0;
/// Maximum rate per minute of necrosis for tumor cells in case of hypoxia with
/// 0 oxygen
constexpr real_t kMaximumNecrosisRate = 1.0 / (6.0 * 60.0);
/// Default oxygen consumption rate of tumor cell
constexpr real_t kDefaultOxygenConsumptionTumorCell = 10.0;
/// Volume parameters
/// Default total volume of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNewTumorCell = 2494.0;
/// Default volume of the nucleus of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNucleusTumorCell = 540.0;
/// Default fraction of fluid volume in a new tumor cell
constexpr real_t kDefaultFractionFluidTumorCell = 0.75;
/// Average time for transformation Random Rate in hours
constexpr real_t kAverageTimeTransformationRandomRate = 38.6;
/// Standard Deviation for transformation Random Rate in hours
constexpr real_t kStandardDeviationTransformationRandomRate = 3.7;
/// Average adhesion time in minutes for Tumor Cell under CAR-T attack before
/// escaping
constexpr real_t kAdhesionTime = 60.0;
/// Min oncoprotein level to be killed by a CAR-T cell
constexpr real_t kOncoproteinLimit = 0.5;
/// Max oncoprotein level
constexpr real_t kOncoproteinSaturation = 2.0;
/// Do not modify this line: difference between saturation and limit
constexpr real_t kOncoproteinDifference =
    kOncoproteinSaturation - kOncoproteinLimit;

/// volume relaxation rate (min^-1) for each state
constexpr real_t kVolumeRelaxationRateAliveCytoplasm = 0.13 / 60.;
constexpr real_t kVolumeRelaxationRateAliveNucleus = 0.22 / 60.;
constexpr real_t kVolumeRelaxationRateAliveFluid = 1.3 / 60.;

constexpr real_t kVolumeRelaxationRateCytoplasmNecroticSwelling = 0.0032 / 60.0;
constexpr real_t kVolumeRelaxationRateNucleusNecroticSwelling = 0.013 / 60.;
constexpr real_t kVolumeRelaxationRateFluidNecroticSwelling = 0.050 / 60.0;

constexpr real_t kVolumeRelaxationRateCytoplasmNecroticLysed = 0.0032 / 60.00;
constexpr real_t kVolumeRelaxationRateNucleusNecroticLysed = 0.013 / 60.;
constexpr real_t kVolumeRelaxationRateFluidNecroticLysed = 0.050 / 60.0;

/// Thresholds in oncoprotein levels for differentiating 4 cancer cell types
constexpr real_t kThresholdCancerCellType1 = 1.5;
constexpr real_t kThresholdCancerCellType2 = 1.0;
constexpr real_t kThresholdCancerCellType3 = 0.5;
constexpr real_t kThresholdCancerCellType4 = 0.0;

///
/// CAR-T Cell Hyperparameters
///

/// Average time in minutes until a CAR-T cell dies
constexpr real_t kAverageMaximumTimeUntillApoptosisCart =
    kDtCycle * 10.0 * 24.0 * 60.0;
/// Default oxygen consumption rate of CAR-T cell
constexpr real_t kDefaultOxygenConsumptionCarT = 1.0;
/// Volume parameters
///  Default total volume of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNewCarTCell = 2494.0;
/// Default volume of the nucleus of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNucleusCarTCell = 540.0;
/// Default fraction of fluid volume in a new CAR-T cell
constexpr real_t kDefaultFractionFluidCarTCell = 0.75;

/// How often a CAR-T cell tries to kill an attached cancer cell in 1/min
constexpr real_t kKillRateCart = 0.06667;
/// How often a CAR-T cell tries to attach to a cancer cell in 1/min
constexpr real_t kAdhesionRateCart = 0.013;
/// Maximum adhesion distance between CAR-T and tumor cells in micrometers
constexpr real_t kMaxAdhesionDistanceCart = 18.0;
/// Minimum adhesion distance between CAR-T and tumor cells in micrometers
constexpr real_t kMinAdhesionDistanceCart = 14.0;

/// Minimum distance in micrometers from the tumor for spawning CAR-T cells
constexpr real_t kMinimumDistanceCarTFromTumor = 50.0;

/// Motility parameters
/// Average persistence time before CAR-T cell moves
constexpr real_t kPersistenceTimeCart = 10;  // 10 minutes
/// Higher bias (\in [0,1]) makes CAR-T movement more directed toward
/// immunostimulatory factor source; while a bias of 0 makes the movement random
constexpr real_t kMigrationBiasCart = 0.5;
/// Migration speed
constexpr real_t kMigrationSpeedCart = 5.0;
/// Elastic constant
constexpr real_t kElasticConstantCart = 0.01;

///
/// Fixed Constants that should not be directly changed
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
/// Number of steps per cycle step Needs to be computed to avoid errors with
/// fmod
constexpr int kStepsPerCycle = kDtCycle / kDt;
/// Steps in one day
constexpr size_t kStepsOneDay = 24 * 60 / kDt;
/// Volume of a single voxel in μm³
constexpr real_t kVoxelVolume =
    (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances) *
    (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances) *
    (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances);
/// 1-kMigrationBiasCart
constexpr real_t kMigrationOneMinusBiasCart = 1.0 - kMigrationBiasCart;
/// Probability of a CAR-T cell to migrate in a given
/// mechanical time step
constexpr real_t kMotilityProbability = kDtMechanics / kPersistenceTimeCart;
/// Probability of a Tumor cell to escape in a given
/// mechanical time step
constexpr real_t kProbabilityEscape = kDtMechanics / (kAdhesionTime + 1e-15);

/// Maximum adhesion distance squared
constexpr real_t kSquaredMaxAdhesionDistanceCart =
    kMaxAdhesionDistanceCart * kMaxAdhesionDistanceCart;
/// Difference between min and max adhesion distance
constexpr real_t kDifferenceCartAdhesionDistances =
    kMaxAdhesionDistanceCart - kMinAdhesionDistanceCart;

/// Radius tumor cell
const real_t gKRadiusTumorCell =
    std::cbrt(kDefaultVolumeNewTumorCell * 3. / (4. * Math::kPi));

/// Radius cart cell
const real_t gKRadiusCarTCell =
    std::cbrt(kDefaultVolumeNewCarTCell * 3. / (4. * Math::kPi));

/// Max Distance for considering two cells as neighbours for force calculations
/// in μm (twice cell radius times kMaxRelativeAdhesionDistance + 0.1 to avoid
/// mismatch because of numerical errors)**2
const real_t gKSquaredMaxDistanceNeighborsForce =
    std::pow(0.1 + 2 * gKRadiusTumorCell * kMaxRelativeAdhesionDistance, 2);

///
/// Some derived values that can also be changed in in the params file
///

/// maximum squared distance to avoid CAR-T pushing
/// tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
/// CAR-T cell will only move to the tumor cell with the adhesion forces
/// (radiusCART + radiusTumorCell + 0.1 to avoid numerical errors)**2
const real_t gKMaxSquaredDistanceCartMovingTowardsTumorCell =
    std::pow(gKRadiusCarTCell + gKRadiusTumorCell + 1, 2);

}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_
