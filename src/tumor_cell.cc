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
#include "tumor_cell.h"

namespace bdm {

TumorCell::TumorCell(const Real3& position) {
  SetPosition(position);
  state_ = TumorCellState::kAlive;  // Default state for new cells
  timer_state_ = 0;  // Initial timer_state

  //volumes
  SetVolume(kDefaultVolumeNewTumorCell);  // Set default volume
  SetFluidFraction(kDefaultFractionFluidTumorCell); // Set default fluid fraction
  SetNuclearVolume(kDefaultVolumeNucleusTumorCell); // Set default nuclear volume
  //target volumes
  SetTargetFractionFluid(kDefaultFractionFluidTumorCell); // Set target fraction of fluid
  SetTargetRelationCytoplasmNucleus((kDefaultVolumeNewTumorCell - kDefaultVolumeNucleusTumorCell) / ( 1e-16 + kDefaultVolumeNucleusTumorCell)); // Set target relation between cytoplasm and nucleus
  SetTargetNucleusSolid(kDefaultVolumeNucleusTumorCell*(1-kDefaultFractionFluidTumorCell)); // Set target nucleus solid volume to real_t
  SetTargetCytoplasmSolid((kDefaultVolumeNewTumorCell - kDefaultVolumeNucleusTumorCell) * (1 - kDefaultFractionFluidTumorCell)); // Set target cytoplasm solid volume to real_t

  SetOncoproteinLevel(SamplePositiveGaussian(kOncoproteinMean,kOncoproteinStandardDeviation)); // Set initial oncoprotein level with a truncated normal distribution 
  SetOncoproteinLevel(0.75);//Debug
  // SetOncoproteinLevel(1.); //Debug
  auto* rm = Simulation::GetActive()->GetResourceManager();
  oxygen_dgrid_ = rm->GetDiffusionGrid("oxygen"); // Pointer to oxygen diffusion grid
  immunostimulatory_factor_dgrid_ = rm->GetDiffusionGrid("immunostimulatory_factor"); // Pointer to immunostimulatory_factor diffusion grid
  SetTransformationRandomRate(); // Set state transition random rate
  attached_to_cart_ = false; // Initially not attached to a cart

  older_velocity_ = {0, 0, 0}; // Initialize the velocity of the cell in the previous step to zero

  //Add Consumption and Secretion
  SetOxygenConsumptionRate(kDefaultOxygenConsumption); // Set default oxygen consumption rate
  SetImmunostimulatoryFactorSecretionRate(kRateSecretionImmunostimulatoryFactor); // Set default immunostimulatory factor secretion rate
  ComputeConstantsConsumptionSecretion(); // Compute constants for all ConsumptionSecretion of Oxygen and Immunostimulatory Factors

}

/// Called when a new agent is created (e.g., after cell division)
void TumorCell::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);

  if (auto* mother = dynamic_cast<TumorCell*>(event.existing_agent)) {//if the cell is created from division
    if (event.GetUid() == CellDivisionEvent::kUid) {
      //Initialize daughter cell from mother cell
      state_ = TumorCellState::kAlive;  // state after division
      timer_state_ = 0;
      //diffusion grids
      oxygen_dgrid_ = mother->oxygen_dgrid_;  // Pointer to the oxygen diffusion grid
      immunostimulatory_factor_dgrid_ = mother->immunostimulatory_factor_dgrid_;  // Pointer to the immunostimulatory_factor diffusion grid
      this->SetOncoproteinLevel(mother->oncoprotein_level_); // inherit oncoprotein level from mother cell
      this->SetOxygenConsumptionRate(mother->GetOxygenConsumptionRate()); // inherit oxygen consumption rate from mother cell
      this->SetImmunostimulatoryFactorSecretionRate(mother->GetImmunostimulatoryFactorSecretionRate()); // inherit immunostimulatory factor secretion rate from mother cell

      // Update the constants for all ConsumptionSecretion
      mother->ComputeConstantsConsumptionSecretion();
      this->ComputeConstantsConsumptionSecretion();


      //divde mother's nuclear volume by 2
      real_t new_nuclear_volume = mother->GetNuclearVolume() / 2.0; // Divide mother's nuclear volume by 2
      mother->SetNuclearVolume(new_nuclear_volume); // Set mother's nuclear volume to the new volume
      this->SetNuclearVolume(new_nuclear_volume);

      //Inherit mother's fluid fraction and velocity
      this->SetFluidFraction(mother->GetFluidFraction()); // Set fluid fraction to mother's fluid fraction
      this->SetOlderVelocity(mother->GetOlderVelocity()); // Copy velocity from mother cell

      //inherit target volumes of the daughter cell
      this->SetTargetFractionFluid(mother->GetTargetFractionFluid());
      this->SetTargetRelationCytoplasmNucleus(mother->GetTargetRelationCytoplasmNucleus());
      this->SetTargetNucleusSolid(mother->GetTargetNucleusSolid());
      this->SetTargetCytoplasmSolid(mother->GetTargetCytoplasmSolid());

      this->SetTransformationRandomRate(); // Set state transition random rate
      this->attached_to_cart_ = false; // Initially not attached to a cart
    }
  }
}

void TumorCell::SetOncoproteinLevel(real_t level) { 
  // std::cout << "Setting oncoprotein level :" << type_ << std::endl;//Debug
  oncoprotein_level_ = level; //oncoprotein_level_
    //cell type
    if (level >= 1.5) {//between 1.5 and 2.0
      type_ = 1;
    } else if (level >= 1.0 && level < 1.5) {
      type_ = 2;
    } else if (level >= 0.5 && level < 1.0) {
      type_ = 3;
    } else {//between 0.0 and 0.5
      type_ = 4;
    }
}

void TumorCell::SetTransformationRandomRate() { 
  // avoid division by zero
  transformation_random_rate_ = 1/(std::max(SamplePositiveGaussian(38.6,3.7)*60., 1e-16));
}

real_t TumorCell::GetTargetTotalVolume() {
  return GetTargetNucleusSolid() * (1 + GetTargetRelationCytoplasmNucleus()) / (1 - GetTargetFractionFluid());
}

// This method explicitly solves the system of exponential relaxation differential equation using a discrete 
// update step. It is used to grow or shrink the volume (and proportions) smoothly toward a desired target 
// volume over time. The relaxations rate controls the speed of convergence
void TumorCell::ChangeVolumeExponentialRelaxationEquation(real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus, real_t relaxation_rate_fluid) {
  // Exponential relaxation towards the target volume
  real_t current_total_volume = GetVolume();
  real_t fluid_fraction= GetFluidFraction();
  real_t nuclear_volume = GetNuclearVolume();

  real_t current_nuclear_solid = nuclear_volume * (1 - fluid_fraction);
  real_t current_cytoplasm_solid = (current_total_volume - nuclear_volume) * (1-fluid_fraction);

  //     std::cout << "time=" << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << 
  // ", current_total_volume=" << current_total_volume <<
  // ", current_nuclear_volume=" << nuclear_volume <<    
  // ", current_cytoplasm_solid=" << current_cytoplasm_solid <<
  // std::endl;

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

//Debug Debug Output params
// std::ofstream file("output/volumes.csv", std::ios::app);
// if (file.is_open()) {

// 	// Write data to CSV file
// 	file << Simulation::GetActive()->GetScheduler()->GetSimulatedSteps()*kDt << ",cytoplasm,"
// 	<< new_volume-total_nuclear << ",nuclear,"
// 	<< total_nuclear <<",fraction fluid,"
// 	<< new_fraction_fluid<< "cytoplasm solid: "
//   <<cytoplasm_solid<<", nuclear solid: "
//   <<nuclear_solid<<"target cytoplasm solid: "
//   <<target_cytoplasm_solid<<", target nuclear solid: "
//   <<target_nucleus_solid_<<", target fluid fraction: "
//   <<target_fraction_fluid_<<", target relation cytoplasm nucleus: "
//   <<target_relation_cytoplasm_nucleus_<< "\n";
// }
// End Debug Output
  
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
Real3 TumorCell::CalculateDisplacement(const InteractionForce* force,
                            real_t squared_radius, real_t dt) {  
  //Debug TumorCell positions
  // Real3 position = GetPosition();
	// std::ofstream file("output/positions_tumor_cell.csv", std::ios::app);
	// if (file.is_open()) {
	// file  << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() << "," 
	// 	<< position[0] << "," << position[1] << "," << position[2] << ", norm:" << position.Norm() << "\n";
	// }
  //End Debug
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a cleaner code instead of hardcoding it here
  squared_radius=kSquaredMaxDistanceNeighborsForce;
  
  Real3 translation_velocity_on_point_mass{0, 0, 0};

  // We check for every neighbor object if they touch us, i.e. push us
  // away and agreagate the velocities

  uint64_t non_zero_neighbor_forces = 0;
  if (!IsStatic()) {
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto calculate_neighbor_forces =
        L2F([&](Agent* neighbor, real_t squared_distance) {
          auto neighbor_force = force->Calculate(this, neighbor);
          if (neighbor_force[0] != 0 || neighbor_force[1] != 0 ||
              neighbor_force[2] != 0) {
            non_zero_neighbor_forces++;
            translation_velocity_on_point_mass[0] += neighbor_force[0];
            translation_velocity_on_point_mass[1] += neighbor_force[1];
            translation_velocity_on_point_mass[2] += neighbor_force[2];
          }
        });
    ctxt->ForEachNeighbor(calculate_neighbor_forces, *this, squared_radius);

    if (non_zero_neighbor_forces > 1) {
      SetStaticnessNextTimestep(false);
    }
  }
  // Two step Adams-Bashforth approximation of the time derivative for position
  // position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 * velocity(t - dt) ]
  movement_at_next_step += translation_velocity_on_point_mass * kDnew + older_velocity_ * kDold;

  older_velocity_ = translation_velocity_on_point_mass;

  // Displacement
  return movement_at_next_step;
}

//Compute new oxygen or immunostimulatory factor concentration after consumption/ secretion
real_t TumorCell::ConsumeSecreteSubstance(int substance_id, real_t old_concentration) {
  // constant1_oxygen_ = 0;  // Debug
  // constant2_oxygen_ = 1.3;  // Debug
  real_t res;
  if (substance_id == oxygen_dgrid_->GetContinuumId()) {
    // consuming oxygen
    res= (old_concentration + constant1_oxygen_) / constant2_oxygen_;
  } else if (substance_id == immunostimulatory_factor_dgrid_->GetContinuumId()) {
    // secreting immunostimulatory factor
    res= (old_concentration + constant1_immunostimulatory_factor_) / constant2_immunostimulatory_factor_;
  } else {
    throw std::invalid_argument("Unknown substance id: " + std::to_string(substance_id));
  }
  return res;
}

void TumorCell::ComputeConstantsConsumptionSecretion() {
  // constant1_= dt · (V_k / V_voxel) · S_k · ρ*_k) 
  // constant2_ = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // where:
  // S_k    = secretion rate of cell k
  // U_k    = uptake (consumption) rate of cell k
  // ρ*_k   = saturation (target) density for secretion
  // V_k    = volume of the cell k
  // V_voxel = volume of the voxel containing the cell
  // dt     = simulation time step
  
  real_t new_volume = GetVolume();
  //compute the constants for the differential equation explicit solution: for oxygen and immunostimulatory factor
  //dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt · (V_k / V_voxel) · S_k · ρ*_k) 
  constant1_oxygen_ = 0.;
  // Scale by the volume of the cell in the Voxel and time step
  constant1_immunostimulatory_factor_ = immunostimulatory_factor_secretion_rate_ * kSaturationDensityImmunostimulatoryFactor * kDtSubstances * (new_volume / kVoxelVolume);
  //1 + dt*(cell_volume/voxel_volume)*(quantity_secretion + quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // Scale by the volume of the cell in the Voxel and time step
  constant2_oxygen_ = 1 + kDtSubstances * (new_volume/ kVoxelVolume) * (oxygen_consumption_rate_);
  // Scale by the volume of the cell in the Voxel and time step
  constant2_immunostimulatory_factor_ = 1 + kDtSubstances * (new_volume/ kVoxelVolume) * (immunostimulatory_factor_secretion_rate_);
}

void TumorCell::StartApoptosis(){

  // std::cout << "Starting apoptosis at time: " << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() <<" IsDead():" << IsDead() << std::endl;//Debug
  // If the cell is already dead, do nothing
  if (IsDead()) {return;}

  // The cell Dies
  SetState(TumorCellState::kApoptotic);

  // Reset timer_state
  SetTimerState(0);  
  // Set target volume to 0 (the cell shrinks)
  SetTargetCytoplasmSolid(0.0);
  SetTargetNucleusSolid(0.0);
  SetTargetFractionFluid(0.0);
  SetTargetRelationCytoplasmNucleus(0.0);
  //Reduce oxygen consumption
  SetOxygenConsumptionRate(GetOxygenConsumptionRate()*kReductionConsumptionDeadCells);
  //Stop Immunostimulatory Factor Secretion
  SetImmunostimulatoryFactorSecretionRate(0.0);
  // Update constants for consumption/secretion differential equation solving
  ComputeConstantsConsumptionSecretion(); 
  // Set type to 5 to indicate dead cell
  SetType(5);
  // std::cout << "Cell type set to 5 (dead) at time: " << Simulation::GetActive()->GetScheduler()->GetSimulatedTime() <<" "<< IsDead()<< std::endl;//Debug
}

/// Main behavior executed at each simulation step
void StateControlGrowProliferate::Run(Agent* agent) {

  auto* sim = Simulation::GetActive();
  if(sim->GetScheduler()->GetSimulatedSteps() % kStepsPerCycle != 0){return;}// Run only every kDtCycle minutes, fmod does not work with the type returned by GetSimulatedTime()
  //Debug
  // // Print simulation minute and number of TumorCell agents
  // int num_steps = sim->GetScheduler()->GetSimulatedSteps();
  // int current_minute = 6 * num_steps;
  // size_t num_cells = sim->GetResourceManager()->GetNumAgents();
  // int current_hour = current_minute / 60;
  // int current_day = current_hour / 24;
  // std::cout << "Dia: " << current_day << "  Hora: " << (current_hour % 24)
  //       << "  Minuto: " << (current_minute % 60)
  //       << "  Numero de celulas: " << num_cells << std::endl;
  // // ----------------------------------------
  // // End Debug

  if (auto* cell = dynamic_cast<TumorCell*>(agent)) {
  // std::cout << "control tumor"<<std::endl;//Debug


    if (cell->IsAttachedToCart()) {
      // If the cell is attached to a cart, skip the state control and growth
      return;
    }
  // std::cout << "contr ol g tum"<<std::endl;//Debug

    // Oxygen levels
    Real3 current_position = cell->GetPosition();
    auto* oxygen_dgrid = cell->GetOxygenDiffusionGrid();  // Pointer to the oxygen diffusion grid
    real_t oxygen_level = oxygen_dgrid->GetValue(current_position);
    // oxygen_level = 30.;  // Debug

    // Debug
    // std::cout << oxygen_level << std::endl;

    real_t multiplier;
    real_t final_rate_transition;

    switch (cell->GetState())
    {
      case TumorCellState::kAlive:{//the cell is growing to real_t its size before mitosis
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);  // Increase timer_state to track time in this state (kDtCycle minutes per step)
        
        if (ShouldEnterNecrosis(oxygen_level, cell)) { // Enter necrosis if oxygen level is too low
          return; // Exit the function to prevent further processing
        }

        //volume change
        cell->ChangeVolumeExponentialRelaxationEquation(kVolumeRelaxarionRateAliveCytoplasm,
                                                        kVolumeRelaxarionRateAliveNucleus,
                                                        kVolumeRelaxarionRateAliveFluid); // The cell grows to real_t its size
        //cell state control
        multiplier = 1.0; // Default multiplier for transition cycle
        if (oxygen_level < kOxygenSaturationInProliferation) {//oxygen threshold for considering an effect on the proliferation cycle
          multiplier = (oxygen_level-kOxygenLimitForProliferation)/(kOxygenSaturationInProliferation-kOxygenLimitForProliferation);
        }
        if(oxygen_level < kOxygenLimitForProliferation) {
          multiplier = 0.0; // If oxygen is below the limit, set multiplier to 0
        }
        // double multiplier1 = multiplier; //Debug


        final_rate_transition= cell->GetTransformationRandomRate() * multiplier * cell->GetOncoproteinLevel(); // Calculate the rate of state change based on oxygen level and oncoprotein (min^-1)

        // //Debug
        // int current_time = sim->GetScheduler()->GetSimulatedSteps()* kDt; // Get the current time step in minutes
        // std::ofstream file("output/simulation_data_mine" + std::to_string(current_time/(12*60)) + ".csv", std::ios::app);
        // if (file.is_open()) {
        // file  << oxygen_level << "," 
        //      << cell->GetOncoproteinLevel() << ","
        //      <<cell->GetTransformationRandomRate()<< "," 
        //      << final_rate_transition << "\n";
        // }
        // //End Debug

            //Debug Debug Output params
        // std::ofstream file2("output/params_o2_oncoprotein.csv", std::ios::app);
        // if (file2.is_open()) {

        //   // Write data to CSV file
        //   file2 << currennt_time << ",multiplier1,"
        //   << multiplier1 << ",multiplier2,"
        //   << multiplier2 << ",transition_rate,"
        // << final_rate_transition
        //   <<"\n";
        // }
        // End Debug Output
        //End Debug

        // //volume change
        // cell->ChangeVolumeExponentialRelaxationEquation(kVolumeRelaxarionRateAliveCytoplasm,
        //                                                 kVolumeRelaxarionRateAliveNucleus,
        //                                                 kVolumeRelaxarionRateAliveFluid); // The cell grows to real_t its size
        // //cell state control
        // multiplier = 1.0; // Default multiplier for transition cycle
        // if (oxygen_level < kOxygenSaturationInProliferation) {//oxygen threshold for considering an effect on the proliferation cycle
        //   multiplier = (oxygen_level-kOxygenLimitForProliferation)/(kOxygenSaturationInProliferation-kOxygenLimitForProliferation);
        // }
        // if(oxygen_level < kOxygenLimitForProliferation) {
        //   multiplier = 0.0; // If oxygen is below the limit, set multiplier to 0
        // }

        
        // final_rate_transition= cell->GetTransformationRandomRate() * multiplier * cell->GetOncoproteinLevel(); // Calculate the rate of state change based on oxygen level and oncoprotein (min^-1)
        
        real_t time_to_wait=1e100; // Set a very large time to avoid division by zero
        if (final_rate_transition > 0) {
          time_to_wait = 1./final_rate_transition; // Calculate the time to transition (in minutes )
        }
        if (time_to_wait< cell->GetTimerState()) { // If the timer_state exceeds the time to transition, change state (this is a fixed duration transition)
          //mitosis: cell divides
          cell->SetState(TumorCellState::kAlive);
          cell->Divide();
          cell->SetTimerState(0);  // Reset timer_state
        }
        break;
      }
      case TumorCellState::kNecroticSwelling:{//the cell is swelling before lysing
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);  // Increase timer_state to track time in this state (kDtCycle minutes per step)
        //volume change
        // The cell swells
        cell->ChangeVolumeExponentialRelaxationEquation(kVolumeRelaxarionRateCytoplasmNecroticSwelling,
                                                        kVolumeRelaxarionRateNucleusNecroticSwelling,
                                                        kVolumeRelaxarionRateFluidNecroticSwelling);
        if (cell->GetVolume() >= 2*kDefaultVolumeNewTumorCell) { // If the cell has swollen to 2 times its original volume, it lyses
          cell->SetState(TumorCellState::kNecroticLysed); // Change state to necrotic lysed
          cell->SetTimerState(0);  // Reset timer_state
          // Set target volume to 0 (the cell will shrink)
          cell->SetTargetCytoplasmSolid(0.0); 
          cell->SetTargetNucleusSolid(0.0); 
          cell->SetTargetFractionFluid(0.0); 
          cell->SetTargetRelationCytoplasmNucleus(0.0);
          // Stop secretion and consumption rate
          // Stop consumption
          cell->SetOxygenConsumptionRate(0.0);
          // Stop secretion
          cell->SetImmunostimulatoryFactorSecretionRate(0.0);
          // Update constants for all ConsumptionSecretion of Oxygen and Immunostimulatory Factors
          cell->ComputeConstantsConsumptionSecretion();
        }
        break;
      }
      case TumorCellState::kNecroticLysed:{//the cell is shirinking and will be removed after a certain time
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);  // Increase timer_state to track time in this state (kDtCycle minutes per step)
        //volume change
        // The cell shrinks
        cell->ChangeVolumeExponentialRelaxationEquation(kVolumeRelaxarionRateCytoplasmNecroticLysed,
                                                        kVolumeRelaxarionRateNucleusNecroticLysed,
                                                        kVolumeRelaxarionRateFluidNecroticLysed);
        if (kTimeLysis < cell->GetTimerState()) { // If the timer_state exceeds the time to transition (this is a fixed duration transition)
          //remove the cell from the simulation
          auto* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      case TumorCellState::kApoptotic:{

        cell->SetTimerState(cell->GetTimerState() + kDtCycle);  // Increase timer_state to track time in this state (kDtCycle minutes per step)
        
        // The cell shrinks
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
        Log::Error("StateControlGrowProliferate::Run", "Unknown TumorCellState");
        break;
      }
    }
  } else {
    Log::Error("StateControlGrowProliferate::Run", "SimObject is not a TumorCell");
  }
}

// computes the probability of the cell entering necrosis
bool StateControlGrowProliferate::ShouldEnterNecrosis(real_t oxygen_level,TumorCell* cell) const {
  //necrosis probability
  real_t multiplier= 0.0; // Default multiplier for necrosis probability
  
  if (oxygen_level < kOxygenLimitForNecrosis){//oxygen threshold for considering necrosis
    multiplier = (kOxygenLimitForNecrosis-oxygen_level)/(kOxygenLimitForNecrosis-kOxygenLimitForNecrosisMaximum);
  }
  if (oxygen_level < kOxygenLimitForNecrosisMaximum) {//threshold for maximum necrosis probability
    multiplier = 1.0;
  }

  real_t probability_necrosis= kDtCycle //multiply by kDtCycle since each timestep is kDtCycle minutes
  * kMaximumNecrosisRate * multiplier; // Calculate the probability of necrosis based on oxygen level

  auto* sim = Simulation::GetActive();
  auto* random = sim->GetRandom();
  bool enter_necrosis = random->Uniform(0, 1) < probability_necrosis;
  if(enter_necrosis){ // If the random number is less than the probability, enter necrosis
    cell->SetState(TumorCellState::kNecroticSwelling); // If oxygen is too low, enter necrosis
    cell->SetTimerState(0);  // Reset timer_state

    //Stop Secretion and reduce consumption
    // Stop secretion
    cell->SetImmunostimulatoryFactorSecretionRate(0.0);
    // Reduce consumption
    cell->SetOxygenConsumptionRate(cell->GetOxygenConsumptionRate()*kReductionConsumptionDeadCells);
    // Update constants for all ConsumptionSecretion of Oxygen and Immunostimulatory Factors
    cell->ComputeConstantsConsumptionSecretion();

    // The cell will swell getting filled with fluid
    cell->SetTargetCytoplasmSolid(0);
    cell->SetTargetNucleusSolid(0);
    cell->SetTargetFractionFluid(1.0); // Set target fraction of fluid to 1.0
    cell->SetTargetRelationCytoplasmNucleus(0.0);
    cell->SetType(5); // Set type to 5 to indicate dead cell
  }
  return enter_necrosis; // Return whether the cell entered necrosis
}

}  // namespace bdm