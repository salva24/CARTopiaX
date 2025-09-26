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

#include "hyperparams.h"
#include "core/real_t.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

namespace bdm {

// Define the static member kUid for SimParam
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

// Function to read a JSON file and load the parameters
// The parameters that are not found in the file are assigned their default
// value
void SimParam::LoadParams(const std::string& filename) {
  nlohmann::json jfile;
  std::ifstream file(filename);
  // If the json can be opened, try to read it, otherwise leave jfile empty
  if (file.is_open()) {
    try {
      file >> jfile;
    } catch (const std::exception& e) {
      std::cerr << "Error reading JSON: " << e.what() << std::endl;
    }
  }

  // Load parameters from JSON file
  if (jfile.contains("kSeed")) {
    kSeed = jfile["kSeed"].get<int>();
  } else {
    kSeed = 1;
  }

  if (jfile.contains("kOutputPerformanceStatistics")) {
    kOutputPerformanceStatistics =
        jfile["kOutputPerformanceStatistics"].get<bool>();
  } else {
    kOutputPerformanceStatistics = false;
  }

  if (jfile.contains("kTotalMinutesToSimulate")) {
    kTotalMinutesToSimulate = jfile["kTotalMinutesToSimulate"].get<int>();
  } else {
    // Default: simulate 30 days
    kTotalMinutesToSimulate = 30 * 24 * 60;
  }

  if (jfile.contains("kInitialRadiusTumor")) {
    kInitialRadiusTumor = jfile["kInitialRadiusTumor"].get<double>();
  } else {
    // Default: initial radius of the tumor is 150 micrometers
    kInitialRadiusTumor = 150.0;
  }

  if (jfile.contains("kBoundedSpaceLength")) {
    kBoundedSpaceLength = jfile["kBoundedSpaceLength"].get<double>();
  } else {
    kBoundedSpaceLength = 1000.0;
  }

  if (jfile.contains("kTreatment")) {
    kTreatment.clear();
    // Parse the JSON object to handle string keys
    for (auto& [key, value] : jfile["kTreatment"].items()) {
      int day = std::stoi(key);
      int dose = value.get<int>();
      kTreatment[day] = dose;
    }
  } else {
    // Example: On day 0 and day 8, 3957 CAR-T cells are introduced (matching
    // the initial tumor cell count). avoid
    // cppcoreguidelines-avoid-magic-numbers warning
    kTreatment = {{0, 3957}, {8, 3957}};
  }

  if (jfile.contains("kDtSubstances")) {
    kDtSubstances = jfile["kDtSubstances"].get<double>();
  } else {
    kDtSubstances = 0.01;
  }

  if (jfile.contains("kDtMechanics")) {
    kDtMechanics = jfile["kDtMechanics"].get<double>();
  } else {
    kDtMechanics = 0.1;
  }

  if (jfile.contains("kDtCycle")) {
    kDtCycle = jfile["kDtCycle"].get<double>();
  } else {
    kDtCycle = 6.0;
  }

  if (jfile.contains("kDt")) {
    kDt = jfile["kDt"].get<double>();
  } else {
    kDt = kDtMechanics;
  }

  if (jfile.contains("kOutputCsvInterval")) {
    kOutputCsvInterval = jfile["kOutputCsvInterval"].get<int>();
  } else {
    kOutputCsvInterval = 12 * 60 / kDt;
  }

  if (jfile.contains("kVolumeRelaxationRateCytoplasmApoptotic")) {
    kVolumeRelaxationRateCytoplasmApoptotic =
        jfile["kVolumeRelaxationRateCytoplasmApoptotic"].get<double>();
  } else {
    kVolumeRelaxationRateCytoplasmApoptotic = 1.0 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateNucleusApoptotic")) {
    kVolumeRelaxationRateNucleusApoptotic =
        jfile["kVolumeRelaxationRateNucleusApoptotic"].get<double>();
  } else {
    kVolumeRelaxationRateNucleusApoptotic = 0.35 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateFluidApoptotic")) {
    kVolumeRelaxationRateFluidApoptotic =
        jfile["kVolumeRelaxationRateFluidApoptotic"].get<double>();
  } else {
    kVolumeRelaxationRateFluidApoptotic = 0.0;
  }

  if (jfile.contains("kTimeApoptosis")) {
    kTimeApoptosis = jfile["kTimeApoptosis"].get<double>();
  } else {
    kTimeApoptosis = 8.6 * 60;
  }

  if (jfile.contains("kReductionConsumptionDeadCells")) {
    kReductionConsumptionDeadCells =
        jfile["kReductionConsumptionDeadCells"].get<double>();
  } else {
    kReductionConsumptionDeadCells = 0.1;
  }

  if (jfile.contains("kResolutionGridSubstances")) {
    kResolutionGridSubstances = jfile["kResolutionGridSubstances"].get<int>();
  } else {
    kResolutionGridSubstances = 50;
  }

  if (jfile.contains("kDiffusionCoefficientOxygen")) {
    kDiffusionCoefficientOxygen =
        jfile["kDiffusionCoefficientOxygen"].get<double>();
  } else {
    kDiffusionCoefficientOxygen = 100000;
  }

  if (jfile.contains("kDecayConstantOxygen")) {
    kDecayConstantOxygen = jfile["kDecayConstantOxygen"].get<double>();
  } else {
    kDecayConstantOxygen = 0.1;
  }

  if (jfile.contains("kDiffusionCoefficientImmunostimulatoryFactor")) {
    kDiffusionCoefficientImmunostimulatoryFactor =
        jfile["kDiffusionCoefficientImmunostimulatoryFactor"].get<double>();
  } else {
    kDiffusionCoefficientImmunostimulatoryFactor = 1000;
  }

  if (jfile.contains("kDecayConstantImmunostimulatoryFactor")) {
    kDecayConstantImmunostimulatoryFactor =
        jfile["kDecayConstantImmunostimulatoryFactor"].get<double>();
  } else {
    kDecayConstantImmunostimulatoryFactor = 0.016;
  }

  if (jfile.contains("kOxygenReferenceLevel")) {
    kOxygenReferenceLevel = jfile["kOxygenReferenceLevel"].get<double>();
  } else {
    kOxygenReferenceLevel = 38.0;
  }

  if (jfile.contains("kInitialOxygenLevel")) {
    kInitialOxygenLevel = jfile["kInitialOxygenLevel"].get<double>();
  } else {
    kInitialOxygenLevel = 38.0;
  }

  if (jfile.contains("kOxygenSaturation")) {
    kOxygenSaturation = jfile["kOxygenSaturation"].get<double>();
  } else {
    kOxygenSaturation = 30.0;
  }

  if (jfile.contains("kRepulsionTumorTumor")) {
    kRepulsionTumorTumor = jfile["kRepulsionTumorTumor"].get<double>();
  } else {
    kRepulsionTumorTumor = 10.0;
  }

  if (jfile.contains("kRepulsionCartCart")) {
    kRepulsionCartCart = jfile["kRepulsionCartCart"].get<double>();
  } else {
    kRepulsionCartCart = 50.0;
  }

  if (jfile.contains("kRepulsionCartTumor")) {
    kRepulsionCartTumor = jfile["kRepulsionCartTumor"].get<double>();
  } else {
    kRepulsionCartTumor = 50.0;
  }

  if (jfile.contains("kRepulsionTumorCart")) {
    kRepulsionTumorCart = jfile["kRepulsionTumorCart"].get<double>();
  } else {
    kRepulsionTumorCart = 10.0;
  }

  if (jfile.contains("kMaxRelativeAdhesionDistance")) {
    kMaxRelativeAdhesionDistance =
        jfile["kMaxRelativeAdhesionDistance"].get<double>();
  } else {
    kMaxRelativeAdhesionDistance = 1.25;
  }

  if (jfile.contains("kAdhesionTumorTumor")) {
    kAdhesionTumorTumor = jfile["kAdhesionTumorTumor"].get<double>();
  } else {
    kAdhesionTumorTumor = 0.4;
  }

  if (jfile.contains("kAdhesionCartCart")) {
    kAdhesionCartCart = jfile["kAdhesionCartCart"].get<double>();
  } else {
    kAdhesionCartCart = 0.0;
  }

  if (jfile.contains("kAdhesionCartTumor")) {
    kAdhesionCartTumor = jfile["kAdhesionCartTumor"].get<double>();
  } else {
    kAdhesionCartTumor = 0.0;
  }

  if (jfile.contains("kAdhesionTumorCart")) {
    kAdhesionTumorCart = jfile["kAdhesionTumorCart"].get<double>();
  } else {
    kAdhesionTumorCart = 0.0;
  }

  if (jfile.contains("kLengthBoxMechanics")) {
    kLengthBoxMechanics = jfile["kLengthBoxMechanics"].get<double>();
  } else {
    kLengthBoxMechanics = 22;
  }

  if (jfile.contains("kDnew")) {
    kDnew = jfile["kDnew"].get<double>();
  } else {
    kDnew = 1.5 * kDtMechanics;
  }

  if (jfile.contains("kDold")) {
    kDold = jfile["kDold"].get<double>();
  } else {
    kDold = -0.5 * kDtMechanics;
  }

  if (jfile.contains("kRateSecretionImmunostimulatoryFactor")) {
    kRateSecretionImmunostimulatoryFactor =
        jfile["kRateSecretionImmunostimulatoryFactor"].get<double>();
  } else {
    kRateSecretionImmunostimulatoryFactor = 10.0;
  }

  if (jfile.contains("kSaturationDensityImmunostimulatoryFactor")) {
    kSaturationDensityImmunostimulatoryFactor =
        jfile["kSaturationDensityImmunostimulatoryFactor"].get<double>();
  } else {
    kSaturationDensityImmunostimulatoryFactor = 1.0;
  }

  if (jfile.contains("kOncoproteinMean")) {
    kOncoproteinMean = jfile["kOncoproteinMean"].get<double>();
  } else {
    kOncoproteinMean = 1.0;
  }

  if (jfile.contains("kOncoproteinStandardDeviation")) {
    kOncoproteinStandardDeviation =
        jfile["kOncoproteinStandardDeviation"].get<double>();
  } else {
    kOncoproteinStandardDeviation = 0.25;
  }

  if (jfile.contains("kOxygenSaturationInProliferation")) {
    kOxygenSaturationInProliferation =
        jfile["kOxygenSaturationInProliferation"].get<double>();
  } else {
    kOxygenSaturationInProliferation = 38.0;
  }

  if (jfile.contains("kOxygenLimitForProliferation")) {
    kOxygenLimitForProliferation =
        jfile["kOxygenLimitForProliferation"].get<double>();
  } else {
    kOxygenLimitForProliferation = 10.0;
  }

  if (jfile.contains("kOxygenLimitForNecrosis")) {
    kOxygenLimitForNecrosis = jfile["kOxygenLimitForNecrosis"].get<double>();
  } else {
    kOxygenLimitForNecrosis = 5.0;
  }

  if (jfile.contains("kOxygenLimitForNecrosisMaximum")) {
    kOxygenLimitForNecrosisMaximum =
        jfile["kOxygenLimitForNecrosisMaximum"].get<double>();
  } else {
    kOxygenLimitForNecrosisMaximum = 2.5;
  }

  if (jfile.contains("kTimeLysis")) {
    kTimeLysis = jfile["kTimeLysis"].get<double>();
  } else {
    kTimeLysis = 60 * 24 * 60.0;
  }

  if (jfile.contains("kMaximumNecrosisRate")) {
    kMaximumNecrosisRate = jfile["kMaximumNecrosisRate"].get<double>();
  } else {
    kMaximumNecrosisRate = 1.0 / (6.0 * 60.0);
  }

  if (jfile.contains("kDefaultOxygenConsumptionTumorCell")) {
    kDefaultOxygenConsumptionTumorCell =
        jfile["kDefaultOxygenConsumptionTumorCell"].get<double>();
  } else {
    kDefaultOxygenConsumptionTumorCell = 10.0;
  }

  if (jfile.contains("kDefaultVolumeNewTumorCell")) {
    kDefaultVolumeNewTumorCell =
        jfile["kDefaultVolumeNewTumorCell"].get<double>();
  } else {
    kDefaultVolumeNewTumorCell = 2494.0;
  }

  if (jfile.contains("kDefaultVolumeNucleusTumorCell")) {
    kDefaultVolumeNucleusTumorCell =
        jfile["kDefaultVolumeNucleusTumorCell"].get<double>();
  } else {
    kDefaultVolumeNucleusTumorCell = 540.0;
  }

  if (jfile.contains("kDefaultFractionFluidTumorCell")) {
    kDefaultFractionFluidTumorCell =
        jfile["kDefaultFractionFluidTumorCell"].get<double>();
  } else {
    kDefaultFractionFluidTumorCell = 0.75;
  }

  if (jfile.contains("kAverageTimeTransformationRandomRate")) {
    kAverageTimeTransformationRandomRate =
        jfile["kAverageTimeTransformationRandomRate"].get<double>();
  } else {
    kAverageTimeTransformationRandomRate = 38.6;
  }

  if (jfile.contains("kStandardDeviationTransformationRandomRate")) {
    kStandardDeviationTransformationRandomRate =
        jfile["kStandardDeviationTransformationRandomRate"].get<double>();
  } else {
    kStandardDeviationTransformationRandomRate = 3.7;
  }

  if (jfile.contains("kAdhesionTime")) {
    kAdhesionTime = jfile["kAdhesionTime"].get<double>();
  } else {
    kAdhesionTime = 60.0;
  }

  if (jfile.contains("kOncoproteinLimit")) {
    kOncoproteinLimit = jfile["kOncoproteinLimit"].get<double>();
  } else {
    kOncoproteinLimit = 0.5;
  }

  if (jfile.contains("kOncoproteinSaturation")) {
    kOncoproteinSaturation = jfile["kOncoproteinSaturation"].get<double>();
  } else {
    kOncoproteinSaturation = 2.0;
  }

  // Difference between saturation and limit. This is always calculated here
  kOncoproteinDifference = kOncoproteinSaturation - kOncoproteinLimit;

  if (jfile.contains("kVolumeRelaxationRateAliveCytoplasm")) {
    kVolumeRelaxationRateAliveCytoplasm =
        jfile["kVolumeRelaxationRateAliveCytoplasm"].get<double>();
  } else {
    kVolumeRelaxationRateAliveCytoplasm = 0.13 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateAliveNucleus")) {
    kVolumeRelaxationRateAliveNucleus =
        jfile["kVolumeRelaxationRateAliveNucleus"].get<double>();
  } else {
    kVolumeRelaxationRateAliveNucleus = 0.22 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateAliveFluid")) {
    kVolumeRelaxationRateAliveFluid =
        jfile["kVolumeRelaxationRateAliveFluid"].get<double>();
  } else {
    kVolumeRelaxationRateAliveFluid = 1.3 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateCytoplasmNecroticSwelling")) {
    kVolumeRelaxationRateCytoplasmNecroticSwelling =
        jfile["kVolumeRelaxationRateCytoplasmNecroticSwelling"].get<double>();
  } else {
    kVolumeRelaxationRateCytoplasmNecroticSwelling = 0.0032 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateNucleusNecroticSwelling")) {
    kVolumeRelaxationRateNucleusNecroticSwelling =
        jfile["kVolumeRelaxationRateNucleusNecroticSwelling"].get<double>();
  } else {
    kVolumeRelaxationRateNucleusNecroticSwelling = 0.013 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateFluidNecroticSwelling")) {
    kVolumeRelaxationRateFluidNecroticSwelling =
        jfile["kVolumeRelaxationRateFluidNecroticSwelling"].get<double>();
  } else {
    kVolumeRelaxationRateFluidNecroticSwelling = 0.050 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateCytoplasmNecroticLysed")) {
    kVolumeRelaxationRateCytoplasmNecroticLysed =
        jfile["kVolumeRelaxationRateCytoplasmNecroticLysed"].get<double>();
  } else {
    kVolumeRelaxationRateCytoplasmNecroticLysed = 0.0032 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateNucleusNecroticLysed")) {
    kVolumeRelaxationRateNucleusNecroticLysed =
        jfile["kVolumeRelaxationRateNucleusNecroticLysed"].get<double>();
  } else {
    kVolumeRelaxationRateNucleusNecroticLysed = 0.013 / 60.0;
  }

  if (jfile.contains("kVolumeRelaxationRateFluidNecroticLysed")) {
    kVolumeRelaxationRateFluidNecroticLysed =
        jfile["kVolumeRelaxationRateFluidNecroticLysed"].get<double>();
  } else {
    kVolumeRelaxationRateFluidNecroticLysed = 0.050 / 60.0;
  }

  if (jfile.contains("kThresholdCancerCellType1")) {
    kThresholdCancerCellType1 =
        jfile["kThresholdCancerCellType1"].get<double>();
  } else {
    kThresholdCancerCellType1 = 1.5;
  }

  if (jfile.contains("kThresholdCancerCellType2")) {
    kThresholdCancerCellType2 =
        jfile["kThresholdCancerCellType2"].get<double>();
  } else {
    kThresholdCancerCellType2 = 1.0;
  }

  if (jfile.contains("kThresholdCancerCellType3")) {
    kThresholdCancerCellType3 =
        jfile["kThresholdCancerCellType3"].get<double>();
  } else {
    kThresholdCancerCellType3 = 0.5;
  }

  if (jfile.contains("kThresholdCancerCellType4")) {
    kThresholdCancerCellType4 =
        jfile["kThresholdCancerCellType4"].get<double>();
  } else {
    kThresholdCancerCellType4 = 0.0;
  }

  if (jfile.contains("kAverageMaximumTimeUntillApoptosisCart")) {
    kAverageMaximumTimeUntillApoptosisCart =
        jfile["kAverageMaximumTimeUntillApoptosisCart"].get<double>();
  } else {
    kAverageMaximumTimeUntillApoptosisCart = kDtCycle * 10.0 * 24.0 * 60.0;
  }

  if (jfile.contains("kDefaultOxygenConsumptionCarT")) {
    kDefaultOxygenConsumptionCarT =
        jfile["kDefaultOxygenConsumptionCarT"].get<double>();
  } else {
    kDefaultOxygenConsumptionCarT = 1.0;
  }

  if (jfile.contains("kDefaultVolumeNewCarTCell")) {
    kDefaultVolumeNewCarTCell =
        jfile["kDefaultVolumeNewCarTCell"].get<double>();
  } else {
    kDefaultVolumeNewCarTCell = 2494.0;
  }

  if (jfile.contains("kKillRateCart")) {
    kKillRateCart = jfile["kKillRateCart"].get<double>();
  } else {
    kKillRateCart = 0.06667;
  }

  if (jfile.contains("kAdhesionRateCart")) {
    kAdhesionRateCart = jfile["kAdhesionRateCart"].get<double>();
  } else {
    kAdhesionRateCart = 0.013;
  }

  if (jfile.contains("kMaxAdhesionDistanceCart")) {
    kMaxAdhesionDistanceCart = jfile["kMaxAdhesionDistanceCart"].get<double>();
  } else {
    kMaxAdhesionDistanceCart = 18.0;
  }

  if (jfile.contains("kMinAdhesionDistanceCart")) {
    kMinAdhesionDistanceCart = jfile["kMinAdhesionDistanceCart"].get<double>();
  } else {
    kMinAdhesionDistanceCart = 14.0;
  }

  if (jfile.contains("kMinimumDistanceCarTFromTumor")) {
    kMinimumDistanceCarTFromTumor =
        jfile["kMinimumDistanceCarTFromTumor"].get<double>();
  } else {
    kMinimumDistanceCarTFromTumor = 50.0;
  }

  if (jfile.contains("kPersistenceTimeCart")) {
    kPersistenceTimeCart = jfile["kPersistenceTimeCart"].get<double>();
  } else {
    kPersistenceTimeCart = 10;
  }

  if (jfile.contains("kMigrationBiasCart")) {
    kMigrationBiasCart = jfile["kMigrationBiasCart"].get<double>();
  } else {
    kMigrationBiasCart = 0.5;
  }

  if (jfile.contains("kMigrationSpeedCart")) {
    kMigrationSpeedCart = jfile["kMigrationSpeedCart"].get<double>();
  } else {
    kMigrationSpeedCart = 5.0;
  }

  if (jfile.contains("kElasticConstantCart")) {
    kElasticConstantCart = jfile["kElasticConstantCart"].get<double>();
  } else {
    kElasticConstantCart = 0.01;
  }

  //
  // Computed constants that should not be directly changed
  //
  // Calculate steps per cycle. This is always calculated here
  kStepsPerCycle = kDtCycle / kDt;
  // Calculate steps per day. This is always calculated here
  kStepsOneDay = 24 * 60 / kDt;
  // Calculate the volume of a single mechanical voxel in μm³
  kVoxelVolume =
      (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances) *
      (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances) *
      (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances);
  // 1-kMigrationBiasCart
  kMigrationOneMinusBiasCart = 1.0 - kMigrationBiasCart;
  // Probability of a CAR-T cell to migrate in a given
  // mechanical time step
  kMotilityProbability = kDtMechanics / kPersistenceTimeCart;
  // Probability of a Tumor cell to escape in a given
  // mechanical time step
  kProbabilityEscape = kDtMechanics / (kAdhesionTime + 1e-15);
  // Maximum adhesion distance squared
  kSquaredMaxAdhesionDistanceCart =
      kMaxAdhesionDistanceCart * kMaxAdhesionDistanceCart;
  // Difference between min and max adhesion distance
  kDifferenceCartAdhesionDistances =
      kMaxAdhesionDistanceCart - kMinAdhesionDistanceCart;
  // Radius tumor cell
  kRadiusTumorCell =
      std::cbrt(kDefaultVolumeNewTumorCell * 3. / (4. * Math::kPi));
  // Radius cart cell
  kRadiusCarTCell =
      std::cbrt(kDefaultVolumeNewCarTCell * 3. / (4. * Math::kPi));
  // Max Distance for considering two cells as neighbours for force calculations
  // in μm (twice cell radius times kMaxRelativeAdhesionDistance + 0.1 to avoid
  // mismatch because of numerical errors)**2
  kSquaredMaxDistanceNeighborsForce =
      std::pow(0.1 + 2 * kRadiusTumorCell * kMaxRelativeAdhesionDistance, 2);

  //
  // Last constant that is derived from others but can be changed directly
  //
  // maximum squared distance to avoid CAR-T pushing
  // tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
  // CAR-T cell will only move to the tumor cell with the adhesion forces
  // (radiusCART + radiusTumorCell + 1 to avoid numerical errors)**2
  kMaxSquaredDistanceCartMovingTowardsTumorCell =
      std::pow(kRadiusCarTCell + kRadiusTumorCell + 1, 2);
}

void SimParam::PrintParams() const {
  std::cout << "\n========================================\n";
  std::cout << "        SIMULATION PARAMETERS\n";
  std::cout << "========================================\n\n";

  ///
  /// General Hyperparameters
  ///
  std::cout << "///\n";
  std::cout << "/// General Hyperparameters\n";
  std::cout << "///\n\n";

  std::cout << "Seed for random number generation: " << kSeed << "\n";
  std::cout << "Output Performance Statistics: "
            << (kOutputPerformanceStatistics ? "true" : "false") << "\n";
  std::cout << "Total simulation time in minutes (30 days): "
            << kTotalMinutesToSimulate << "\n";
  std::cout << "Initial radius of the spherical tumor (micrometers): "
            << kInitialRadiusTumor << "\n";
  std::cout << "Length of the bounded space (micrometers): "
            << kBoundedSpaceLength << "\n\n";

  /// Treatment Dosages
  std::cout << "/// Treatment Dosages\n";
  std::cout << "///\n\n";
  std::cout << "CAR-T cell infusion schedule:\n";
  if (kTreatment.empty()) {
    std::cout << "  No treatment scheduled\n";
  } else {
    for (const auto& [day, dose] : kTreatment) {
      std::cout << "  Day " << day << ": " << dose << " CAR-T cells\n";
    }
  }
  std::cout << "\n";

  /// Time steps
  std::cout << "/// Time steps\n";
  std::cout << "///\n\n";
  std::cout << "Time step for substances secretion/consumption (minutes): "
            << kDtSubstances << "\n";
  std::cout << "Time step for the cell mechanics (minutes): " << kDtMechanics
            << "\n";
  std::cout << "Time step for the cell cycle (minutes): " << kDtCycle << "\n";
  std::cout << "General time step for the simulation: " << kDt << "\n";
  std::cout << "Output CSV interval: " << kOutputCsvInterval << "\n\n";

  /// Apoptotic cells volume change
  std::cout << "/// Apoptotic cells volume change\n";
  std::cout << "///\n\n";
  std::cout << "Volume relaxation rate cytoplasm apoptotic: "
            << kVolumeRelaxationRateCytoplasmApoptotic << "\n";
  std::cout << "Volume relaxation rate nucleus apoptotic: "
            << kVolumeRelaxationRateNucleusApoptotic << "\n";
  std::cout << "Volume relaxation rate fluid apoptotic: "
            << kVolumeRelaxationRateFluidApoptotic << "\n";
  std::cout << "Time in minutes until an apoptotic cell is removed: "
            << kTimeApoptosis << "\n";
  std::cout << "Reduction of consumption rate of dead cells when they enter "
               "necrosis: "
            << kReductionConsumptionDeadCells << "\n\n";

  /// Chemicals
  ///
  std::cout << "/// Chemicals\n";
  std::cout << "///\n\n";
  std::cout << "Number of voxels per axis for the substances grid: "
            << kResolutionGridSubstances << "\n";
  std::cout << "Diffusion coefficient of oxygen (μm²/min): "
            << kDiffusionCoefficientOxygen << "\n";
  std::cout << "Decay constant of oxygen (min⁻¹): " << kDecayConstantOxygen
            << "\n";
  std::cout << "Diffusion coefficient of immunostimulatory factor (μm²/min): "
            << kDiffusionCoefficientImmunostimulatoryFactor << "\n";
  std::cout << "Decay constant of immunostimulatory factor (min⁻¹): "
            << kDecayConstantImmunostimulatoryFactor << "\n";
  std::cout << "Reference level of oxygen at the boundaries (mmHg): "
            << kOxygenReferenceLevel << "\n";
  std::cout << "Initial oxygen concentration in each voxel (mmHg): "
            << kInitialOxygenLevel << "\n";
  std::cout << "Oxygen saturation in the microenvironment (mmHg): "
            << kOxygenSaturation << "\n\n";

  /// Forces
  ///
  std::cout << "/// Forces\n";
  std::cout << "///\n\n";
  std::cout << "Repulsion coefficient between tumor cells: "
            << kRepulsionTumorTumor << "\n";
  std::cout << "Repulsion coefficient between CAR-T cells: "
            << kRepulsionCartCart << "\n";
  std::cout << "Repulsion coefficient between CAR-T cells and tumor cells: "
            << kRepulsionCartTumor << "\n";
  std::cout << "Repulsion coefficient between tumor cells and CAR-T cells: "
            << kRepulsionTumorCart << "\n";
  std::cout << "Maximum relative adhesion distance for cell interactions: "
            << kMaxRelativeAdhesionDistance << "\n";
  std::cout << "Adhesion coefficient between tumor cells: "
            << kAdhesionTumorTumor << "\n";
  std::cout << "Adhesion coefficient between CAR-T cells: " << kAdhesionCartCart
            << "\n";
  std::cout << "Adhesion coefficient between CAR-T cells and tumor cells: "
            << kAdhesionCartTumor << "\n";
  std::cout << "Adhesion coefficient between tumor cells and CAR-T cells: "
            << kAdhesionTumorCart << "\n";
  std::cout << "Box length for mechanics calculations (micrometers): "
            << kLengthBoxMechanics << "\n";
  std::cout
      << "Coefficient for the current time step in Adams-Bashforth method: "
      << kDnew << "\n";
  std::cout
      << "Coefficient for the previous time step in Adams-Bashforth method: "
      << kDold << "\n";
  std::cout << "Maximum squared distance to avoid CAR-T pushing tumor cells: "
            << kMaxSquaredDistanceCartMovingTowardsTumorCell << "\n\n";

  ///
  /// TumorCell Hyperparameters
  ///
  std::cout << "///\n";
  std::cout << "/// TumorCell Hyperparameters\n";
  std::cout << "///\n\n";

  std::cout << "Rate of secretion of immunostimulatory factor of tumor cells "
               "per minute: "
            << kRateSecretionImmunostimulatoryFactor << "\n";
  std::cout
      << "Saturation density of immunostimulatory factor for tumor cells: "
      << kSaturationDensityImmunostimulatoryFactor << "\n";
  std::cout << "Mean level of oncoprotein expression in tumor cells: "
            << kOncoproteinMean << "\n";
  std::cout << "Standard deviation of oncoprotein expression in tumor cells: "
            << kOncoproteinStandardDeviation << "\n";
  std::cout << "Oxygen saturation level in tumor cells for proliferation: "
            << kOxygenSaturationInProliferation << "\n";
  std::cout << "Limit of oxygen level for tumor cell proliferation: "
            << kOxygenLimitForProliferation << "\n";
  std::cout << "Limit of oxygen to start causing necrosis: "
            << kOxygenLimitForNecrosis << "\n";
  std::cout << "Limit of oxygen to maximum necrosis probability: "
            << kOxygenLimitForNecrosisMaximum << "\n";
  std::cout << "Time in minutes until a lysed necrotic cell is removed: "
            << kTimeLysis << "\n";
  std::cout << "Maximum rate per minute of necrosis for tumor cells: "
            << kMaximumNecrosisRate << "\n";
  std::cout << "Default oxygen consumption rate of tumor cell: "
            << kDefaultOxygenConsumptionTumorCell << "\n";

  std::cout << "\nVolume parameters:\n";
  std::cout << "Default total volume of a new tumor cell (μm³): "
            << kDefaultVolumeNewTumorCell << "\n";
  std::cout << "Default volume of the nucleus of a new tumor cell (μm³): "
            << kDefaultVolumeNucleusTumorCell << "\n";
  std::cout << "Default fraction of fluid volume in a new tumor cell: "
            << kDefaultFractionFluidTumorCell << "\n";

  std::cout << "\nTransformation and adhesion parameters:\n";
  std::cout << "Average time for transformation Random Rate (hours): "
            << kAverageTimeTransformationRandomRate << "\n";
  std::cout << "Standard Deviation for transformation Random Rate (hours): "
            << kStandardDeviationTransformationRandomRate << "\n";
  std::cout
      << "Average adhesion time for Tumor Cell under CAR-T attack (minutes): "
      << kAdhesionTime << "\n";
  std::cout << "Min oncoprotein level to be killed by a CAR-T cell: "
            << kOncoproteinLimit << "\n";
  std::cout << "Max oncoprotein level: " << kOncoproteinSaturation << "\n";
  std::cout << "Oncoprotein difference (saturation - limit): "
            << kOncoproteinDifference << "\n";

  std::cout << "\nVolume relaxation rates for alive cells:\n";
  std::cout << "Volume relaxation rate alive cytoplasm (min⁻¹): "
            << kVolumeRelaxationRateAliveCytoplasm << "\n";
  std::cout << "Volume relaxation rate alive nucleus (min⁻¹): "
            << kVolumeRelaxationRateAliveNucleus << "\n";
  std::cout << "Volume relaxation rate alive fluid (min⁻¹): "
            << kVolumeRelaxationRateAliveFluid << "\n";

  std::cout << "\nVolume relaxation rates for necrotic swelling:\n";
  std::cout << "Volume relaxation rate cytoplasm necrotic swelling (min⁻¹): "
            << kVolumeRelaxationRateCytoplasmNecroticSwelling << "\n";
  std::cout << "Volume relaxation rate nucleus necrotic swelling (min⁻¹): "
            << kVolumeRelaxationRateNucleusNecroticSwelling << "\n";
  std::cout << "Volume relaxation rate fluid necrotic swelling (min⁻¹): "
            << kVolumeRelaxationRateFluidNecroticSwelling << "\n";

  std::cout << "\nVolume relaxation rates for necrotic lysed:\n";
  std::cout << "Volume relaxation rate cytoplasm necrotic lysed (min⁻¹): "
            << kVolumeRelaxationRateCytoplasmNecroticLysed << "\n";
  std::cout << "Volume relaxation rate nucleus necrotic lysed (min⁻¹): "
            << kVolumeRelaxationRateNucleusNecroticLysed << "\n";
  std::cout << "Volume relaxation rate fluid necrotic lysed (min⁻¹): "
            << kVolumeRelaxationRateFluidNecroticLysed << "\n";

  std::cout << "\nThresholds in oncoprotein levels for differentiating cancer "
               "cell types:\n";
  std::cout << "Threshold Cancer Cell Type 1: " << kThresholdCancerCellType1
            << "\n";
  std::cout << "Threshold Cancer Cell Type 2: " << kThresholdCancerCellType2
            << "\n";
  std::cout << "Threshold Cancer Cell Type 3: " << kThresholdCancerCellType3
            << "\n";
  std::cout << "Threshold Cancer Cell Type 4: " << kThresholdCancerCellType4
            << "\n\n";

  ///
  /// CAR-T Cell Hyperparameters
  ///
  std::cout << "///\n";
  std::cout << "/// CAR-T Cell Hyperparameters\n";
  std::cout << "///\n\n";

  std::cout << "Average time in minutes until a CAR-T cell dies: "
            << kAverageMaximumTimeUntillApoptosisCart << "\n";
  std::cout << "Default oxygen consumption rate of CAR-T cell: "
            << kDefaultOxygenConsumptionCarT << "\n";

  std::cout << "\nVolume parameters:\n";
  std::cout << "Default total volume of a new CAR-T cell (μm³): "
            << kDefaultVolumeNewCarTCell << "\n";

  std::cout << "\nKilling and adhesion rates:\n";
  std::cout << "How often a CAR-T cell tries to kill an attached cancer cell "
               "(1/min): "
            << kKillRateCart << "\n";
  std::cout
      << "How often a CAR-T cell tries to attach to a cancer cell (1/min): "
      << kAdhesionRateCart << "\n";
  std::cout << "Maximum adhesion distance between CAR-T and tumor cells "
               "(micrometers): "
            << kMaxAdhesionDistanceCart << "\n";
  std::cout << "Minimum adhesion distance between CAR-T and tumor cells "
               "(micrometers): "
            << kMinAdhesionDistanceCart << "\n";

  std::cout << "\nSpawning parameters:\n";
  std::cout << "Minimum distance from the tumor for spawning CAR-T cells "
               "(micrometers): "
            << kMinimumDistanceCarTFromTumor << "\n";

  std::cout << "\nMotility parameters:\n";
  std::cout << "Average persistence time before CAR-T cell moves: "
            << kPersistenceTimeCart << "\n";
  std::cout << "Migration bias (higher values = more directed movement): "
            << kMigrationBiasCart << "\n";
  std::cout << "Migration speed: " << kMigrationSpeedCart << "\n";
  std::cout << "Elastic constant: " << kElasticConstantCart << "\n\n";

  std::cout << "========================================\n";
  std::cout << "      END OF PARAMETERS SUMMARY\n";
  std::cout << "========================================\n\n";
}

}  // namespace bdm