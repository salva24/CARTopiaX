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

#include "cart_cell.h"
#include "core/environment/uniform_grid_environment.h"

namespace bdm {

CartCell::CartCell(const Real3& position) {
  SetPosition(position);
  // Default state for new cells
  state_ = CartCellState::kAlive;
  // Initial timer_state for apoptotic state
  timer_state_ = 0;

  //volumes
  // Set default volume
  SetVolume(kDefaultVolumeNewCartCell);
  // Set default fluid fraction
  SetFluidFraction(kDefaultFractionFluidCartCell);
  // Set default nuclear volume
  SetNuclearVolume(kDefaultVolumeNucleusCartCell);


  // Pointer to oxygen diffusion grid
  auto* rm = Simulation::GetActive()->GetResourceManager();
  oxygen_dgrid_ = rm->GetDiffusionGrid("oxygen");
  // Pointer to immunostimulatory_factor diffusion grid
  immunostimulatory_factor_dgrid_ = rm->GetDiffusionGrid("immunostimulatory_factor");
  // Initially not attached to a tumor cell
  attached_to_tumor_cell_ = false; 
  // Initialize attached cell pointer
  attached_cell_ptr_ = nullptr;

  // Initialize the velocity of the cell in the previous step to zero
  older_velocity_ = {0, 0, 0};


  SetCurrentLiveTime(kAverageMaximumTimeUntillApoptosisCart);

  //Add Consumption and Secretion
  // Set default oxygen consumption rate
  SetOxygenConsumptionRate(kDefaultOxygenConsumptionCarT);
  // Compute constants for all ConsumptionSecretion of Oxygen
  ComputeConstantsConsumptionSecretion();

}

// Cart cells can move if they are alive and not attached to a tumor cell
bool CartCell::DoesCellMove() {
  return state_ == CartCellState::kAlive && !attached_to_tumor_cell_; 
}


real_t CartCell::GetTargetTotalVolume() {
  return GetTargetNucleusSolid() * (1 + GetTargetRelationCytoplasmNucleus()) / (1 - GetTargetFractionFluid());
}

// This method explicitly solves the system of exponential relaxation differential equation using a discrete 
// update step. It is used to shrink the volume (and proportions) smoothly toward a desired target
// volume over time whe the cell is apoptotic. The relaxations rate controls the speed of convergence
void CartCell::ChangeVolumeExponentialRelaxationEquation(real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus, real_t relaxation_rate_fluid) {
  // Exponential relaxation towards the target volume
  real_t current_total_volume = GetVolume();
  real_t fluid_fraction= GetFluidFraction();
  real_t nuclear_volume = GetNuclearVolume();

  real_t current_nuclear_solid = nuclear_volume * (1 - fluid_fraction);
  real_t current_cytoplasm_solid = (current_total_volume - nuclear_volume) * (1-fluid_fraction);

  real_t current_fluid = fluid_fraction * current_total_volume;

  // Update fluid volume
  real_t new_fluid = current_fluid + kDtCycle* relaxation_rate_fluid * (GetTargetFractionFluid() * current_total_volume - current_fluid);
  // Clamp to zero to prevent negative volumes
  if (new_fluid < 0.0) { new_fluid = 0.0; }
  
  real_t nuclear_fluid = new_fluid* ( nuclear_volume/ current_total_volume);
  // real_t cytoplasm_fluid = new_fluid - nuclear_fluid;

  real_t nuclear_solid = current_nuclear_solid + kDtCycle * relaxation_rate_nucleus * (GetTargetNucleusSolid() - current_nuclear_solid);
  // Clamp to zero to prevent negative volumes
  if (nuclear_solid < 0.0) { nuclear_solid = 0.0; }

  real_t target_cytoplasm_solid = GetTargetRelationCytoplasmNucleus() * GetTargetNucleusSolid();
  real_t cytoplasm_solid = current_cytoplasm_solid + kDtCycle * relaxation_rate_cytoplasm * (target_cytoplasm_solid - current_cytoplasm_solid);
  // Clamp to zero to prevent negative volumes
  if (cytoplasm_solid < 0.0) { cytoplasm_solid = 0.0; }

  real_t new_total_solid= nuclear_solid + cytoplasm_solid;

  real_t total_nuclear= nuclear_solid + nuclear_fluid;

  // real_t total_cytoplasm= cytoplasm_solid + cytoplasm_fluid;

  real_t new_volume = new_total_solid + new_fluid;

  // Avoid division by zero
  real_t new_fraction_fluid = new_fluid / (1e-16 + new_volume);
  
  // Update the cell's properties
  // if the volume has changed
  if (new_volume!= current_total_volume){
    SetVolume(new_volume);
    // Update constants for all ConsumptionSecretion of Oxygen and Immunostimulatory Factors
    ComputeConstantsConsumptionSecretion();
  }
  SetFluidFraction(new_fraction_fluid);
  SetNuclearVolume(total_nuclear);
}

//compute Displacement
Real3 CartCell::CalculateDisplacement(const InteractionForce* force,
                            real_t squared_radius, real_t dt) {

  //Debug CAR-T positions
  // Real3 position = GetPosition();
	// std::ofstream file("output/positions_cart.csv", std::ios::app);
	// if (file.is_open()) {
	// file  << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << "," 
	// 	<< position[0] << "," << position[1] << "," << position[2] << ", norm:" << position.Norm() << "\n";
	// }
  //End Debug


  auto* sim = Simulation::GetActive();

  // real_t h = dt;
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a cleaner code instead of hardcoding it here
  squared_radius=kSquaredMaxDistanceNeighborsForce;

  // the physics force to move the point mass

  Real3 translation_velocity_on_point_mass{0, 0, 0};

  //--------------------------------------------
  //CAR-T self motility (in case of migration)
  //--------------------------------------------
  Real3 current_position = GetPosition();
  auto* ctxt = sim->GetExecutionContext();
  auto* rng= sim->GetRandom();
  Real3 motility;
  if (DoesCellMove()){
    //compute motility
    // if (true) {
    if (rng->Uniform(0.0, 1.0) < kMotilityProbability) {//Debug uncomment
      //random direction as unitary vector
      Real3 random_direction = GenerateRandomDirection();
      Real3 direction_to_immunostimulatory_factor;
      // returns normalized gradient towards the immunostimulatory factor source
      immunostimulatory_factor_dgrid_->GetGradient(current_position, &direction_to_immunostimulatory_factor,true);  
      // motility = bias * direction_to_immunostimulatory_factor + (1-bias)*random_direction
      motility = kMigrationBiasCart * direction_to_immunostimulatory_factor + kMigrationOneMinusBiasCart * random_direction;
      // Convert to unit direction
      if (motility[0]*motility[0] + motility[1]*motility[1] + motility[2]*motility[2] > 0)
        motility.Normalize();
      // Scale by migration speed and add to the velocity
      translation_velocity_on_point_mass += motility * kMigrationSpeedCart;//Debug uncomment
      // std::cout << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() <<"velocity afet motility: " << translation_velocity_on_point_mass[0] <<", "<< translation_velocity_on_point_mass[1] <<", "<< translation_velocity_on_point_mass[2] <<", norm "<< translation_velocity_on_point_mass.Norm() << std::endl;//Debug
    }
    

  }

  //--------------------------------------------
  //CAR-T killing or victim cell escaping
  //--------------------------------------------
  if (state_ == CartCellState::kAlive) {//If cell is not apoptotic

    if (attached_to_tumor_cell_) {//attached to tumor cell
      // try to kill the cancer cell and in case of failure see if it manages to scape
      // the order needs to be this one because it should try to kill before seeing if it scapes 
      if (TryToInduceApoptosis(attached_cell_ptr_,rng) || rng->Uniform(0.0, 1.0) < kProbabilityEscape) {//Debug Uncomment 
        //the cancer cell is detached
        attached_cell_ptr_->SetAttachedToCart(false);
        // std::cout << "Detaching at time: " << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << " IsDead():" << attached_cell->IsDead() << std::endl;//Debug
        //empty ID
        attached_cell_ptr_ = nullptr;
        attached_to_tumor_cell_ = false;
      }
    }

    //--------------------------------------------
    //CAR-T adhesion to victim cell
    //--------------------------------------------
    //Compute forces between the cells and check for a new attachment
    auto calculate_forces_and_elastic_displacement =
      L2F([&](Agent* neighbor, real_t squared_distance) {
        //Adhesion repulsion forces between cells
        // We check for every neighbor object if they touch us, i.e. push us
        // away and aggregate the velocities
        auto neighbor_force = force->Calculate(this, neighbor);
        // std::cout << "Neighbor force: " << neighbor_force[0] << std::endl;//Debug
        translation_velocity_on_point_mass[0] += neighbor_force[0];
        translation_velocity_on_point_mass[1] += neighbor_force[1];
        translation_velocity_on_point_mass[2] += neighbor_force[2];

        //CAR-T adhesion to new victim cell
        Real3 displac = neighbor->GetPosition()-current_position;

        if (TumorCell* cancer_cell = dynamic_cast<TumorCell*>(neighbor)) {
          // if(cancer_cell->GetBoxIdx() == this->GetBoxIdx()){//Debug//debug
          // std::cout <<Simulation::GetActive()->GetScheduler()->GetSimulatedTime() <<"Movement towards tumor cell: " << displac[0]*kElasticConstantCart << std::endl;//Debug
          //movement towards the tumor cells
          real_t sq_norm_displac = displac[0]*displac[0] + displac[1]*displac[1] + displac[2]*displac[2];
          
          //The cart moves towards the tumor cell only if they are not touching already
          //If they are too close the only force affecting is the adhesion force to avoid 
          //CAR-T non-stop pushing tumor cells. In case of being closer than kMaxSquaredDistanceCartMovingTowardsTumorCell
          //there is a probability kProbabilityPushing for the CAR-T to keep pushing the tumor cell
          if (sq_norm_displac > kMaxSquaredDistanceCartMovingTowardsTumorCell) {// || rng->Uniform(0.0, 1.0) < kProbabilityPushing) {//Debug
            translation_velocity_on_point_mass[0] += displac[0] * kElasticConstantCart;
            translation_velocity_on_point_mass[1] += displac[1] * kElasticConstantCart;
            translation_velocity_on_point_mass[2] += displac[2] * kElasticConstantCart;
          }

          //If the CAR-T has not succeeded in attaching to a tumor cell yet, it tries again
          if (!attached_to_tumor_cell_)//Debug uncomment
            TryToGetAttachedTo(cancer_cell, sq_norm_displac, rng);
        // }//Debug
        }
      });
    ctxt->ForEachNeighbor(calculate_forces_and_elastic_displacement, *this, squared_radius);
    // std::cout << "Translation velocity on point mass: " << squared_radius << std::endl;//Debug
  }

  //--------------------------------------------
  // Two step Adams-Bashforth approximation of the time derivative for position
  // position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 * velocity(t - dt) ]
  //--------------------------------------------
  movement_at_next_step += translation_velocity_on_point_mass * kDnew + older_velocity_ * kDold;


  older_velocity_ = translation_velocity_on_point_mass;

  // Displacement
  return movement_at_next_step;
}

// Try to get attached to a tumor cell
void CartCell::TryToGetAttachedTo(TumorCell* victim, real_t squared_distance, Random* rng){
    // std::cout << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << std::endl;//Debug
  // std::cout << "Trying to attach at time: " << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << " is dead: " << victim->IsDead()<< " type: " << victim->GetType() << std::endl;//Debug
  // //Debug
  // real_t current_time = Simulation::GetActive()->GetScheduler()->GetSimulatedTime(); // Get the current time step in minutes
  // std::ofstream file("output/cart_adhesion.csv", std::ios::app);
  // if (file.is_open()) {
  // file  << current_time<< " probability "<<  <<"\n";
  // }
  // //End Debug

  //If the tumor cell is not already attached to a CAR-T cell, is not dead and is not too far away.
  if(!victim->IsAttachedToCart()&& !victim->IsDead() && squared_distance < kSquaredMaxAdhesionDistanceCart) {

    //factor of how high is the oncoprotein level of the cancer cell
    real_t oncoprotein_scale_factor = (victim->GetOncoproteinLevel()-kOncoproteinLimit)/kOncoproteinDifference;
    // Clamp scale_factor to be in [0,1]
    if( oncoprotein_scale_factor > 1.0 )
      oncoprotein_scale_factor = 1.0;
    // If oncoprotein level is lower than the limit the cancer cell does not get detected
    if( oncoprotein_scale_factor <= 0.0 ){
      oncoprotein_scale_factor = 0.0;
      // the probability is going to be 0 so return the function is the most efficient
      return;
    }

    //factor of how far the cancer cell is
    real_t distance_scale_factor= (kMaxAdhesionDistanceCart-std::sqrt(squared_distance))/kDifferenceCartAdhesionDistances;
    //Clamp scale_factor to be in [0,1]. We already checked that it is > 0 because squared_distance < kSquaredMaxAdhesionDistanceCart 
    if( distance_scale_factor > 1.0 )
      distance_scale_factor = 1.0;

    
    // distance_scale_factor=1;//Debug
      
    // //Debug
    // real_t current_time = Simulation::GetActive()->GetScheduler()->GetSimulatedTime(); // Get the current time step in minutes
    // std::ofstream file("output/cart_adhesion.csv", std::ios::app);
    // if (file.is_open()) {
    // file  << current_time<< " probability "<< kAdhesionRateCart * oncoprotein_scale_factor * distance_scale_factor * kDtMechanics <<"\n";
    // }
    //End Debug

    // //Debug
    // #pragma omp critical
    // {
    //   acumulator_probabilities += kAdhesionRateCart * oncoprotein_scale_factor * distance_scale_factor * kDtMechanics;
    // }
    // // ENd Debug


    // It tries to attach the CAR-T cell to the tumor cell with probability kAdhesionRateCart * oncoprotein_scale_factor * distance_scale_factor * kDtMechanics
    if (rng->Uniform(0.0, 1.0) < kAdhesionRateCart * oncoprotein_scale_factor * distance_scale_factor * kDtMechanics) {//Debug uncomment
      //avoid race condition. Only one cell can be attached to the tumor cell.
      #pragma omp critical
      {
        //We need to check again if the victim is not attached to a CAR-T cell yet. This could be made more efficiently with a semaphore for each cancer cell
        if (!victim->IsAttachedToCart()){
        // std::cout << "Attaching at time: " << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << std::endl;//Debug
          attached_to_tumor_cell_ = true;
          attached_cell_ptr_ = victim->GetAgentPtr<TumorCell>();
          victim->SetAttachedToCart(true);
        }
      }
    }
  }
    
  return;
}

//Try to induce apoptosis
bool CartCell::TryToInduceApoptosis(bdm::AgentPointer<TumorCell> attached_cell, Random* rng) {
  // If there is no attached cell (this should not happen)

  if (!attached_to_tumor_cell_)
    return false;

  //factor of how high is the oncoprotein levelof the cancer cell
	real_t scale_factor = (attached_cell->GetOncoproteinLevel()-kOncoproteinLimit)/kOncoproteinDifference;
	// Clamp scale_factor to be in [0,1]
  if( scale_factor > 1.0 )
	  scale_factor = 1.0;
  // If oncoprotein level is lower than the limit the cancer cell does not become apoptotic
  if( scale_factor < 0.0 ){
	  scale_factor = 0.0;
    // the probability is going to be 0 so return the function is the most efficient
    return false;
  }
  //CAR-T success of killing probability: aggressive cancer cells (high oncoprotein level) are more likely to be killed
  bool succeeded =  rng->Uniform(0.0, 1.0) < kKillRateCart * scale_factor * kDtMechanics;
  
    //Debug
    // real_t current_time = Simulation::GetActive()->GetScheduler()->GetSimulatedTime(); // Get the current time step in minutes
    // std::ofstream file("output/cart_killing.csv", std::ios::app);
    // if (file.is_open()) {
    // file  << current_time<< " probability "<< kKillRateCart * scale_factor * kDtMechanics <<"\n";
    // }
    // //End Debug

  //The CAR-T has succeeded to induce apoptosis on the Cancer Cell
	if(succeeded)
    attached_cell->StartApoptosis();

  // std::cout << "Cancer cell apoptosis started at time: " << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << " IsDead():" << attached_cell->IsDead() << std::endl;

  return succeeded;
}

//Compute new oxygen or immunostimulatory factor concentration after consumption/ secretion
real_t CartCell::ConsumeSecreteSubstance(int substance_id, real_t old_concentration) {
  real_t res;
  if (substance_id == oxygen_dgrid_->GetContinuumId()) {
    // consuming oxygen
    res= (old_concentration + constant1_oxygen_) / constant2_oxygen_;
  } else if (substance_id == immunostimulatory_factor_dgrid_->GetContinuumId()) {
    //CAR-T do not change immunostimulatory factor levels
    res= old_concentration;

  } else {
    throw std::invalid_argument("Unknown substance id: " + std::to_string(substance_id));
  }
  return res;
}

//Recompute Consumption constants whenever oxygen_consumption_rate_ or the volume changes
void CartCell::ComputeConstantsConsumptionSecretion() {

  // constant1_= dt · (V_k / V_voxel) · S_k · ρ*_k) 
  // constant2_ = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // where:
  // S_k    = secretion rate of cell k
  // U_k    = uptake (consumption) rate of cell k
  // ρ*_k   = saturation (target) density for secretion
  // V_k    = volume of the cell k
  // V_voxel = volume of the voxel containing the cell
  // dt     = simulation time step
  real_t volume = GetVolume();
  //compute the constants for the differential equation explicit solution: for oxygen and immunostimulatory factor
  //dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt · (V_k / V_voxel) · S_k · ρ*_k) 
  constant1_oxygen_ = 0.;
  //1 + dt*(cell_volume/voxel_volume)*(quantity_secretion + quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // Scale by the volume of the cell in the Voxel and time step
  constant2_oxygen_ = 1 + kDtSubstances * (volume/kVoxelVolume) * (oxygen_consumption_rate_);
}

/// Main behavior executed at each simulation step
void StateControlCart::Run(Agent* agent) {

  auto* sim = Simulation::GetActive();
  if(sim->GetScheduler()->GetSimulatedSteps() % kStepsPerCycle != 0){return;}// Run only every kDtCycle minutes, fmod does not work with the type returned by GetSimulatedTime()

  if (auto* cell = dynamic_cast<CartCell*>(agent)) {

    switch (cell->GetState())
    {
      case CartCellState::kAlive:{//the cell is growing to real_t its size before mitosis

        if (sim->GetRandom()->Uniform(1.0) < kDtCycle/std::max(cell->GetCurrentLiveTime(), 1e-10)) { // Probability of death= 1/CurrentLiveTime, avoiding division by 0
          //the cell Dies
          cell->SetState(CartCellState::kApoptotic);
          // Reset timer_state, it should be 0 anyway
          cell->SetTimerState(0);  
          // Set target volume to 0 (the cell will shrink)
          cell->SetTargetCytoplasmSolid(0.0); 
          cell->SetTargetNucleusSolid(0.0); 
          cell->SetTargetFractionFluid(0.0); 
          cell->SetTargetRelationCytoplasmNucleus(0.0);
          //Reduce oxygen consumption
          cell->SetOxygenConsumptionRate(cell->GetOxygenConsumptionRate()*kReductionConsumptionDeadCells);
          // Update constants for all Consumption of oxygen
          cell->ComputeConstantsConsumptionSecretion(); 
          // Detach from tumor cell if it was attached
          if (cell->IsAttachedToTumorCell()) {
            cell->GetAttachedCellPointer()->SetAttachedToCart(false);
            cell->SetAttachedCellPointer(nullptr);
            cell->SetAttachedToTumorCell(false);
          }
        } else{
          // decrease current life time
          cell->SetCurrentLiveTime((cell->GetCurrentLiveTime() - (kDtCycle*kDtCycle)));
        }
        break;
      }
      case CartCellState::kApoptotic:{
        cell->SetTimerState(cell->GetTimerState() + kDtCycle); 

        cell->ChangeVolumeExponentialRelaxationEquation(kVolumeRelaxarionRateCytoplasmApoptotic,
                                                kVolumeRelaxarionRateNucleusApoptotic,
                                                kVolumeRelaxarionRateFluidApoptotic);
        if (kTimeApoptosis < cell->GetTimerState()) { // If the timer_state exceeds the time to transition (this is a fixed duration transition)
          //remove the cell from the simulation
          auto* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      default:{
        Log::Error("StateControlCart::Run", "Unknown CartCellState");
        break;
      }
    }
  } else {
    Log::Error("StateControlCart::Run", "SimObject is not a CartCell");
  }
}

}  // namespace bdm