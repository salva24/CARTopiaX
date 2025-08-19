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
#include "utils_aux.h"
#include "cart_cell.h"


namespace bdm {

// Samples a Gaussian value with given mean and standard deviation but all negative values are mapped to zero
real_t SamplePositiveGaussian(float mean, float sigma) {
  auto* random = Simulation::GetActive()->GetRandom();
  real_t value = random->Gaus(mean, sigma);
  if(value < 0.) {value = 0.;}
  return value;
}


std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius) {
    // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
    real_t cell_radius = std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi)/2;

    std::vector<Real3> positions;

    // Hexagonal close-packing spacing
    real_t spacing_x = cell_radius * std::sqrt(3.0);
    real_t spacing_y = cell_radius * 2.0;
    real_t spacing_z = cell_radius * std::sqrt(3.0);

    int zc = 0;
    for (real_t z = -sphere_radius; z < sphere_radius; z += spacing_z, ++zc) {
        int xc = 0;
        for (real_t x = -sphere_radius; x < sphere_radius; x += spacing_x, ++xc) {
            int yc = 0;
            for (real_t y = -sphere_radius; y < sphere_radius; y += spacing_y, ++yc) {

                // Compute cell center with HCP offset
                real_t px = x + (zc % 2) * 0.5 * cell_radius;
                real_t py = y + (xc % 2) * cell_radius;
                real_t pz = z;

                real_t dist = std::sqrt(px * px + py * py + pz * pz);

                if (dist <= sphere_radius) {
                    positions.push_back({px, py, pz});
                }
            }
        }
    }

    return positions;
}

//Function to compute the number of tumor cells of each type and the radius of the tumor
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, size_t, real_t> ComputeNumberTumorCellsAndRadius() {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  size_t total_num_tumor_cells = 0;
  size_t num_tumor_cells_type1 = 0;
  size_t num_tumor_cells_type2 = 0;
  size_t num_tumor_cells_type3 = 0;
  size_t num_tumor_cells_type4 = 0;
  size_t num_tumor_cells_type5_dead = 0;
  size_t num_alive_cart = 0;

  real_t max_dist_sq = 0.0;

  rm->ForEachAgent([&](const Agent* agent) {
    if (auto* tumor_cell = dynamic_cast<const TumorCell*>(agent)) {
      total_num_tumor_cells++;
      const auto& pos = agent->GetPosition();
      real_t dist_sq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
      if (dist_sq > max_dist_sq) {
        max_dist_sq = dist_sq;
      }

      // Count tumor cells by type
      switch (tumor_cell->GetType()) {
        case 1: num_tumor_cells_type1++; break;
        case 2: num_tumor_cells_type2++; break;
        case 3: num_tumor_cells_type3++; break;
        case 4: num_tumor_cells_type4++; break;
        case 5: num_tumor_cells_type5_dead++; break;
        default: break;
      }
    } else if (auto* cart_cell = dynamic_cast<const CartCell*>(agent)) {
      if (cart_cell->GetState() == CartCellState::kAlive)
        num_alive_cart++;
    }
  });
  return {total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, num_alive_cart,std::sqrt(max_dist_sq)};
}

// Function to spawn CAR-T cell dosages
void SpawnCart::operator()() {
  auto* simulation = Simulation::GetActive();
  auto* scheduler = simulation->GetScheduler();

  //This function only executes each day
  if (scheduler->GetSimulatedSteps() % frequency_ != 0)
    return;
  //See if there is any dosage to apply in this day
  size_t current_day = scheduler->GetSimulatedTime()/ (60.0*24);
  size_t cells_to_spawn=0;

  if (kTreatment.find(current_day) != kTreatment.end())
    cells_to_spawn = kTreatment.at(current_day);

  //if there are cells to spawn in the treatment
  if (cells_to_spawn > 0) {

    //compute tumor radius
    auto* rm = simulation->GetResourceManager();
    real_t max_dist_sq = 0.0;

    rm->ForEachAgent([&](const Agent* agent) {
      if (std::strcmp(agent->GetTypeName(), "TumorCell") == 0) {
        const auto& pos = agent->GetPosition();
        real_t dist_sq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
        if (dist_sq > max_dist_sq) {
          max_dist_sq = dist_sq;
        }
      }
    });

    //the car-t spawns at least 50 micrometers away from the tumor
    real_t minimum_squared_radius= std::sqrt(max_dist_sq)+50;
    minimum_squared_radius *=minimum_squared_radius;

    //for generating car-t positions
    auto* rng= simulation->GetRandom();
    auto* param = simulation->GetParam();
    real_t min_b= param->min_bound;
    real_t max_b= param->max_bound;
    real_t px,py,pz,radi_sq;

    auto* ctxt = simulation->GetExecutionContext();
    CartCell* cart = nullptr;

    for( unsigned int i=0 ;i < cells_to_spawn ; i++ ){
      //look for a valid position
      do{
          px = rng->Uniform(min_b, max_b);
          py = rng->Uniform(min_b, max_b);
          pz = rng->Uniform(min_b, max_b);
          radi_sq = px*px + py*py + pz*pz;
      }while(minimum_squared_radius > radi_sq);

      //spawn CAR-T
      cart = new CartCell({px,py,pz});
      cart->AddBehavior(new StateControlCart());
      ctxt->AddAgent(cart);
    }
	}

}

// Function to output summary CSV
void OutputSummary::operator()() {

  auto* simulation = Simulation::GetActive();
  auto* scheduler = simulation->GetScheduler();
  auto current_step = scheduler->GetSimulatedSteps();

  if (current_step % frequency_ == 0) {
    std::ofstream file("output/final_data.csv", std::ios::app);
    if (file.is_open()) {
      if (current_step == 0) {
        file << "total_days,total_hours,total_minutes,tumor_radius,num_cells,num_tumor_cells,tumor_cells_type1,tumor_cells_type2,tumor_cells_type3,tumor_cells_type4,tumor_cells_type5_dead,num_alive_cart\n";// Header for CSV file
      }

      // Calculate time in days, hours, minutes
      double total_minutes = scheduler->GetSimulatedTime();
      double total_hours = total_minutes / 60.0;
      double total_days = total_hours / 24.0;

      // Count total cells, tumor cells of each type and tumor radius
      size_t total_num_tumor_cells;
      size_t num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead,num_alive_cart;
      real_t tumor_radius;
      std::tie(total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, num_alive_cart, tumor_radius) = ComputeNumberTumorCellsAndRadius();      
      real_t total_num_cells=simulation->GetResourceManager()->GetNumAgents();

      //If a dosage is administred this exact time the numbers are not seen in 
      //the resource manager yet because of how BioDynaMo is built
      //therefore we need to add the just added new cells to the statistics here.
      size_t current_day_int = static_cast<size_t>(total_days);
      if (current_step % kStepsOneDay==0 && kTreatment.find(current_day_int) != kTreatment.end()){
        size_t just_spawned_cells= kTreatment.at(current_day_int);
        total_num_cells += just_spawned_cells;
        num_alive_cart += just_spawned_cells;
      }
      

      // Write data to CSV file
      file << total_days << ","
       << total_hours << ","
       << total_minutes << ","
       << tumor_radius << ","
       << total_num_cells << ","
       << total_num_tumor_cells << ","
       << num_tumor_cells_type1 << ","
       << num_tumor_cells_type2 << ","
       << num_tumor_cells_type3 << ","
       << num_tumor_cells_type4 << ","
       << num_tumor_cells_type5_dead << ","
       << num_alive_cart << "\n";
    }
  }
}

//Function to generate a random direction unitary vector
Real3 GenerateRandomDirection() {
  auto* rnd= Simulation::GetActive()->GetRandom();
  // Generate a random direction vector
  double theta = kTwicePi * rnd->Uniform(0.0, 1.0);
  double phi = Math::kPi * rnd->Uniform(0.0, 1.0);

  double sin_phi = std::sin(phi);
  double cos_phi = std::cos(phi);

  return {sin_phi*std::cos(theta), sin_phi*std::sin(theta), cos_phi};
}

}  // namespace bdm