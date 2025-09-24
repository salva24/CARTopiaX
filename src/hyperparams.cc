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
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

namespace bdm {

// Define the static member kUid for SimParam
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

// Function to read a JSON file and load the parameters
// The parameters that are not found in the file are assigned their default value
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
    }
    std::cout << "Using seed: " << kSeed << std::endl;


























    // // maximum squared distance to avoid CAR-T pushing
    // // tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
    // // CAR-T cell will only move to the tumor cell with the adhesion forces
    // // (radiusCART + radiusTumorCell + 0.1 to avoid numerical errors)**2
    // sparams->kMaxSquaredDistanceCartMovingTowardsTumorCell =
    //     std::pow(sparams->kRadiusCarTCell + sparams->kRadiusTumorCell + 1, 2);






}

}  // namespace bdm