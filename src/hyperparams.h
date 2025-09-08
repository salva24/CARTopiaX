/*
 * Copyright 2025 compiler-research.org, Salvador de la Torre Gonzalez, Luciana Melina Luque
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

#include <map>
#include <string>

#include "core/real_t.h"

namespace bdm {

///This file contains hyperparameters used in the simulation. Change: In a future version it needs to be changed into a params file with no need to be recompiled

///
/// TumorCell Hyperparameters
///

/// Rate of secretion of immunostimulatory factor of tumor cells per minute
constexpr real_t kRateSecretionImmunostimulatoryFactor= 10.0; 
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
constexpr real_t kOxygenLimitForNecrosisMaximum= 2.5; 
/// Time in minutes until a lysed necrotic cell is removed from the simulation
constexpr real_t kTimeLysis = 60*24*60.; 
/// Rate of cell division in min**-1
constexpr real_t kDivisionRate = 0.02717 / 60.0; 
/// Maximum rate per minute of necrosis for tumor cells in case of hypoxia with 0 oxygen
constexpr real_t kMaximumNecrosisRate= 1.0 / (6.0 * 60.0); 
/// Default oxygen consumption rate of tumor cell
constexpr real_t kDefaultOxygenConsumption = 10.0; 
///Volume parameters
/// Default total volume of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNewTumorCell = 2494.0; 
/// Default volume of the nucleus of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNucleusTumorCell = 540.0; 
/// Default fraction of fluid volume in a new tumor cell
constexpr real_t kDefaultFractionFluidTumorCell = 0.75; 
/// Average adhesion time in minutes for Tumor Cell under CAR-T attack before escaping
constexpr real_t kAdhesionTime = 60.0;
/// Min oncoprotein level to be killed by a CAR-T cell
constexpr real_t kOncoproteinLimit = 0.5;
/// Max oncoprotein level
constexpr real_t kOncoproteinSaturation = 2.0;
///Do not modify this line: difference between saturation and limit
constexpr real_t kOncoproteinDifference = kOncoproteinSaturation - kOncoproteinLimit;

///volume relaxation rate (min^-1) for each state
constexpr real_t kVolumeRelaxarionRateAliveCytoplasm =0.13/60.;// 0.27/ 60.0;
constexpr real_t kVolumeRelaxarionRateAliveNucleus = 0.22/60.;//0.33/60.
constexpr real_t kVolumeRelaxarionRateAliveFluid = 1.3/60.;//3.0/60.

constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticSwelling = 0.0032/60.0;
constexpr real_t kVolumeRelaxarionRateNucleusNecroticSwelling = 0.013/60.;
constexpr real_t kVolumeRelaxarionRateFluidNecroticSwelling = 0.050/60.0;

constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticLysed = 0.0032/60.00;
constexpr real_t kVolumeRelaxarionRateNucleusNecroticLysed = 0.013/60.;
constexpr real_t kVolumeRelaxarionRateFluidNecroticLysed = 0.050/60.0;

///
/// General Hyperparameters
///

/// Seed for random number generation
constexpr int kSeed =3; 

/// Output Performance Statistics
constexpr bool kOutputPerformanceStatistics = true;

/// 0.01 minutes time step for substances secretion/consumption
constexpr real_t kDtSubstances = 0.01; 
/// 0.1 minutes time step for the cell mechanics
constexpr real_t kDtMechanics = 0.1; 
/// 6 minutes time step for the cell cycle
constexpr real_t kDtCycle = 6.0; 

/// General time step for the simulation: it is the same as kDtMechanics, do not modify this line
constexpr real_t kDt = kDtMechanics; 
/// Number of steps per cycle step, do not modify this line. Needs to be computed to avoid errors with fmod
constexpr int kStepsPerCycle = kDtCycle / kDt; 

/// Output little summary each half a day
constexpr int kOutputCsvInterval = 12*60/kDt;


/// Total simulation time in minutes (30 days)
constexpr int kTotalMinutesToSimulate = 500; //30*24*60
/// Length of the bounded space in micrometers
constexpr int kBoundedSpaceLength = 1000; 
/// Initial radius of the spherical tumor (group of cancer cells) in micrometers
constexpr real_t kInitialRadiusTumor = 150; 

///Do not modify this line: Twice Pi
constexpr real_t kTwicePi=2.*Math::kPi;


constexpr real_t kVolumeRelaxarionRateCytoplasmApoptotic = 1.0/60.0;
constexpr real_t kVolumeRelaxarionRateNucleusApoptotic = 0.35/60.0;
constexpr real_t kVolumeRelaxarionRateFluidApoptotic = 0.0;
/// Time in minutes until an apoptotic cell is removed from the simulation
constexpr real_t kTimeApoptosis = 8.6*60; 
/// Reduction of consumption rate of dead cells when they enter necrosis
constexpr real_t kReductionConsumptionDeadCells= 0.1; 



///Chemicals
/// Number of voxels per axis for the substances grid
constexpr int kResolutionGridSubstances = 50; //50 // voxels per axis
/// Volume of a single voxel in μm³ (do not modify this line)
constexpr real_t kVoxelVolume = (kBoundedSpaceLength / kResolutionGridSubstances)*(kBoundedSpaceLength / kResolutionGridSubstances)*(kBoundedSpaceLength/ kResolutionGridSubstances); //Do not modify this line
/// Diffusion coefficient of oxygen in μm²/min
constexpr real_t kDiffusionCoefficientOxygen = 100000; // 100000 micrometers^2/minute
/// Decay constant of oxygen in min⁻¹
constexpr real_t kDecayConstantOxygen = 0.1; // 0.1 minutes^-1
/// Diffusion coefficient of immunostimulatory factor in μm²/min
constexpr real_t kDiffusionCoefficientImmunostimulatoryFactor = 1000; // 1000 micrometers^2/minute
/// Decay constant of immunostimulatory factor in min⁻¹
constexpr real_t kDecayConstantImmunostimulatoryFactor = 0.016; // 0.016 minutes^-1
/// Time step for oxygen diffusion in minutes
constexpr real_t kTimeStepOxygen = 0.0005; // 0.001 minutes CHANGE
/// Time step for immunostimulatory factor diffusion in minutes
constexpr real_t kTimeStepImmunostimulatoryFactor = 0.01; // 0.01 minutes
/// Reference level of oxygen at the boundaries in mmHg
constexpr real_t kOxygenReferenceLevel = 38.; // Reference level of oxygen at the boundaries of the simulation space in mmHg
/// Initial oxygen concentration in each voxel in mmHg
constexpr real_t kInitialOxygenLevel = 38.0; // Initial voxel concentration of oxygen in mmHg
/// Oxygen saturation in the microenvironment in mmHg
constexpr real_t kOxygenSaturation = 30.0; //30.0 // Oxygen saturation in mmHg in microenvironment
///Forces
/// Repulsion coeficient between tumor cells
constexpr real_t kRepulsionTumorTumor = 10.0; 
/// Repulsion coeficient between CAR-T cells
constexpr real_t kRepulsionCartCart = 50.0;   
/// Repulsion coeficient between CAR-T cells and tumor cells
constexpr real_t kRepulsionCartTumor = 50.0;  
/// Repulsion coeficient between tumor cells and CAR-T cells
constexpr real_t kRepulsionTumorCart = 10.0;  
/// Maximum relative adhesion distance for cell interactions
constexpr real_t kMaxRelativeAdhesionDistance =1.25; 
/// Adhesion coeficient between tumor cells
constexpr real_t kAdhesionTumorTumor = 0.4; 
/// Adhesion coeficient between CAR-T cells
constexpr real_t kAdhesionCartCart = 0.0;   
/// Adhesion coeficient between CAR-T cells and tumor cells
constexpr real_t kAdhesionCartTumor = 0.0;  
/// Adhesion coeficient between tumor cells and CAR-T cells
constexpr real_t kAdhesionTumorCart = 0.0;  

///Do not change
//coefficientes for the two step Adams-Bashforth approximation of the time derivative for position 
//position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 * velocity(t - dt) ]
/// Coefficient for the current time step in the Adams-Bashforth method (dt * 1.5)
constexpr real_t kDnew = 1.5 * kDtMechanics; 
/// Coefficient for the previous time step in the Adams-Bashforth method (dt * -0.5)
constexpr real_t kDold = -0.5 * kDtMechanics;

///Do not change this line
const real_t kLengthBoxMechanics =22; // Length of the box for mechanics in micrometers

///
/// CAR-T Cell Hyperparameters
///
constexpr real_t kAverageMaximumTimeUntillApoptosisCart= kDtCycle* 10.0 * 24.0 * 60.0;
/// Default oxygen consumption rate of CAR-T cell
constexpr real_t kDefaultOxygenConsumptionCarT = 1.0; 
///Volume parameters
/// Default total volume of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNewCartCell = 2494.0; 
/// Default volume of the nucleus of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNucleusCartCell = 540.0; 
/// Default fraction of fluid volume in a new CAR-T cell
constexpr real_t kDefaultFractionFluidCartCell = 0.75; 

/// How often a CAR-T cell tries to kill an attached cancer cell
constexpr real_t kKillRateCart = 0.06667; // 1/min
/// How often a CAR-T cell tries to attach to a cancer cell
constexpr real_t kAdhesionRateCart = 0.01;//0.2;//0.0092//0.0085; // 1/min
/// Maximum adhesion distance between CAR-T and tumor cells
constexpr real_t kMaxAdhesionDistanceCart = 18.0;//micrometers
/// Minimum adhesion distance between CAR-T and tumor cells
constexpr real_t kMinAdhesionDistanceCart = 14.0;//micrometers

/// Motility parameters
/// Average persistence time before CAR-T cell moves
constexpr real_t kPersistenceTimeCart = 10; // 10 minutes
///Higher bias (\in [0,1]) makes CAR-T movement more directed toward immunostimulatory factor source; while a bias of 0 makes the movement random
constexpr real_t kMigrationBiasCart = 1.;//0.5;//DEBUG
/// Migration speed
constexpr real_t kMigrationSpeedCart = 5.0;
///Elastic constant
constexpr real_t kElasticConstantCart = 0.01;

/// Probability of a CAR-T cell to keep pushing a tumor cell when they are very close
// constexpr real_t kProbabilityPushing = 0.2;//Debug

/// Treatment Dosages
///
/// Specifies the CAR-T cell infusion schedule as a map where:
///   - The key represents the day of treatment (starting from day 0).
///   - The value represents the number of CAR-T cells administered on that day.
/// Example: On day 0 and day 8, 3964 CAR-T cells are introduced (matching the initial tumor cell count).
inline std::map<size_t, size_t> kTreatment = {//Debug
    {0, 0},  // Day 0: administer 3957 CAR-T cells
    {8, 0}   // Day 8: administer 3957 CAR-T cells
};

/// Do not modify this line:  1-kMigrationBiasCart
constexpr real_t kMigrationOneMinusBiasCart = 1.0 - kMigrationBiasCart;
/// Do not modify this line:  probability of a CAR-T cell to migrate in a given mechanical time step
constexpr real_t kMotilityProbability = kDtMechanics / kPersistenceTimeCart;
/// Do not modify this line:  probability of a Tumor cell to escape in a given mechanical time step
constexpr real_t kProbabilityEscape = kDtMechanics / ( kAdhesionTime + 1e-15 );
/// Do not modify this line: Steps in one day
constexpr size_t kStepsOneDay = 24*60/kDt;
/// Do not modify this line: maximum adhesion distance squared
constexpr real_t kSquaredMaxAdhesionDistanceCart = kMaxAdhesionDistanceCart*kMaxAdhesionDistanceCart;
/// Do not modify this line: difference between min and max adhesion distance
constexpr real_t kDifferenceCartAdhesionDistances = kMaxAdhesionDistanceCart - kMinAdhesionDistanceCart;


///Do not change this line: radius tumor cell
const real_t kRadiusTumorCell = std::cbrt(kDefaultVolumeNewTumorCell * 3. / (4. *  Math::kPi));

///Do not change this line: radius cart cell
const real_t kRadiusCartCell = std::cbrt(kDefaultVolumeNewCartCell * 3. / (4. *  Math::kPi));

//Do not change this line: maximum squared distance to avoid CAR-T pushing tumor cells
//If a CAR-T and a Tumor Cell are closer than this distance, the CAR-T cell will only move to the tumor cell with the adhesion forces
//(radiusCART + radiusTumorCell + 0.1 to avoid numerical errors)**2
const real_t kMaxSquaredDistanceCartMovingTowardsTumorCell = std::pow(kRadiusCartCell + kRadiusTumorCell + 1, 2);

///Max Distance for considering two cells as neighbours for force calculations in μm
///Do not change this line
const real_t kSquaredMaxDistanceNeighborsForce = std::pow( 0.1+ 2* std::max(kRadiusTumorCell, kRadiusCartCell) * kMaxRelativeAdhesionDistance,2);// (twice biggest cell radius (in case to cells tha maximum size encounter each other) times kMaxRelativeAdhesionDistance + 0.1 to avoid mismatch because of numerical errors)**2



}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_
