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

#include "agents/tumor_cell.h"
#include "params/hyperparams.h"
#include "utils/utils_aux.h"
#include "core/agent/agent.h"
#include "core/agent/cell_division_event.h"
#include "core/agent/new_agent_event.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/functor.h"
#include "core/interaction_force.h"
#include "core/param/param.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include "core/simulation.h"
#include "core/util/log.h"
#include "core/util/random.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace bdm {

TumorCell::TumorCell(const Real3& position) {
  SetPosition(position);
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  // volumes
  // Set default volume
  SetVolume(SamplePositiveGaussian(sparams->default_volume_new_tumor_cell,sparams->std_volume_new_tumor_cell));
  // Set default fluid fraction
  SetFluidFraction(sparams->default_fraction_fluid_tumor_cell);
  // Set default nuclear volume
  SetNuclearVolume(sparams->default_volume_nucleus_tumor_cell);
  // target volumes
  // Set target fraction of fluid
  SetTargetFractionFluid(sparams->default_fraction_fluid_tumor_cell);
  // Set target relation between cytoplasm and nucleus
  SetTargetRelationCytoplasmNucleus(
      (sparams->default_volume_new_tumor_cell -
       sparams->default_volume_nucleus_tumor_cell) /
      (kEpsilon + sparams->default_volume_nucleus_tumor_cell));
  // Set target nucleus solid volume
  SetTargetNucleusSolid(sparams->default_volume_nucleus_tumor_cell *
                        (1 - sparams->default_fraction_fluid_tumor_cell));
  // Set target cytoplasm solid volume
  SetTargetCytoplasmSolid((sparams->default_volume_new_tumor_cell -
                           sparams->default_volume_nucleus_tumor_cell) *
                          (1 - sparams->default_fraction_fluid_tumor_cell));

  // Set initial oncoprotein level with a truncated normal distribution
  SetOncoproteinLevel(SamplePositiveGaussian(
      sparams->oncoprotein_mean, sparams->oncoprotein_standard_deviation));
  ResourceManager* rm = Simulation::GetActive()->GetResourceManager();
  // Pointer to oxygen diffusion grid
  oxygen_dgrid_ = rm->GetDiffusionGrid("oxygen");
  // Pointer to immunostimulatory_factor diffusion grid if it is added in the simulation
  immunostimulatory_factor_dgrid_ =
    sparams->add_immunostimulatory_factor
        ? rm->GetDiffusionGrid("immunostimulatory_factor")
        : nullptr;
  glucose_dgrid_ = sparams->add_glucose
        ? rm->GetDiffusionGrid("glucose")//Debug
        : nullptr;
  // Set state transition random rate
  SetTransformationRandomRate();
  // Set basal death probability
  SetBasalDeathProbability(sparams->basal_necrosis_probability_cancer_cells);

  // Add Consumption and Secretion
  // Set default oxygen consumption rate
  SetOxygenConsumptionRate(sparams->default_oxygen_consumption_tumor_cell);
  // Set default immunostimulatory factor secretion rate
  SetImmunostimulatoryFactorSecretionRate(
      sparams->rate_secretion_immunostimulatory_factor);
  SetGlucoseConsumptionRate(sparams->default_glucose_consumption_tumor_cell);
  // Compute constants for all ConsumptionSecretion of Oxygen, glucose and
  // Immunostimulatory Factor
  ComputeConstantsConsumptionSecretion();
}

/// Called when a new agent is created (e.g., after cell division)
void TumorCell::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);

  // if the cell is created from division
  if (auto* mother = dynamic_cast<TumorCell*>(event.existing_agent)) {
    if (event.GetUid() == CellDivisionEvent::kUid) {
      // Initialize daughter cell from mother cell
      // state after division
      state_ = TumorCellState::kAlive;
      timer_state_ = 0;
      // diffusion grids
      // Pointer to the oxygen diffusion grid
      oxygen_dgrid_ = mother->oxygen_dgrid_;
      // Pointer to the immunostimulatory_factor diffusion grid
      immunostimulatory_factor_dgrid_ = mother->IsImmunostimulatoryFactorDefined()
                                           ? mother->immunostimulatory_factor_dgrid_
                                           : nullptr;
      // Pointer to the glucose diffusion grid
      glucose_dgrid_ = mother->IsGlucoseDefined() ? mother->glucose_dgrid_ : nullptr;
      // inherit oncoprotein level from mother cell
      this->SetOncoproteinLevel(mother->oncoprotein_level_);
      // inherit oxygen consumption from mother cell
      this->SetOxygenConsumptionRate(mother->GetOxygenConsumptionRate());
      // inherit immunostimulatory factor secretion rate from mother cell
      this->SetImmunostimulatoryFactorSecretionRate(
          mother->GetImmunostimulatoryFactorSecretionRate());
      // inherit glucose consumption from mother cell
      this->SetGlucoseConsumptionRate(mother->GetGlucoseConsumptionRate());

      // Update the constants for all ConsumptionSecretion
      mother->ComputeConstantsConsumptionSecretion();
      this->ComputeConstantsConsumptionSecretion();

      // divide mother's nuclear volume by 2
      const real_t new_nuclear_volume = mother->GetNuclearVolume() / kHalf;
      // Set mother's nuclear volume to the new volume
      mother->SetNuclearVolume(new_nuclear_volume);
      this->SetNuclearVolume(new_nuclear_volume);

      // Inherit mother's fluid fraction and velocity
      // Set fluid fraction to mother's fluid fraction
      this->SetFluidFraction(mother->GetFluidFraction());
      // Copy velocity from mother cell
      this->SetOlderVelocity(mother->GetOlderVelocity());

      // inherit target volumes of the daughter cell
      this->SetTargetFractionFluid(mother->GetTargetFractionFluid());
      this->SetTargetRelationCytoplasmNucleus(
          mother->GetTargetRelationCytoplasmNucleus());
      this->SetTargetNucleusSolid(mother->GetTargetNucleusSolid());
      this->SetTargetCytoplasmSolid(mother->GetTargetCytoplasmSolid());

      // Set state transition random rate
      this->SetTransformationRandomRate();
      // Set basal death probability
      this->SetBasalDeathProbability(mother->GetBasalDeathProbability());
      // Initially not attached to a cart
      this->attached_to_cart_ = false;
    }
  }
}

void TumorCell::SetOncoproteinLevel(real_t level) {
  oncoprotein_level_ = level;
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  // cell type
  if (level >= sparams->threshold_cancer_cell_type1) {
    // between 1.5 and 2.0
    SetType(TumorCellType::kType1);
  } else if (level >= sparams->threshold_cancer_cell_type2 &&
             level < sparams->threshold_cancer_cell_type1) {
    SetType(TumorCellType::kType2);
  } else if (level >= sparams->threshold_cancer_cell_type3 &&
             level < sparams->threshold_cancer_cell_type2) {
    SetType(TumorCellType::kType3);
  } else if (level >= sparams->threshold_cancer_cell_type4 &&
             level < sparams->threshold_cancer_cell_type3) {
    SetType(TumorCellType::kType4);
  } else {
    // undefined type
    SetType(TumorCellType::kType0);
  }
}

void TumorCell::SetTransformationRandomRate() {
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  // avoid division by zero with kEpsilon
  transformation_random_rate_ =
      1 /
      (std::max(SamplePositiveGaussian(
                    sparams->average_time_transformation_random_rate,
                    sparams->standard_deviation_transformation_random_rate) *
                    kMinutesInAnHour,
                kEpsilon));
}

real_t TumorCell::GetTargetTotalVolume() const {
  return GetTargetNucleusSolid() * (1 + GetTargetRelationCytoplasmNucleus()) /
         (1 - GetTargetFractionFluid());
}

// This method explicitly solves the system of exponential relaxation
// differential equation using a discrete update step. It is used to grow or
// shrink the volume (and proportions) smoothly toward a desired target volume
// over time. The relaxations rate controls the speed of convergence
void TumorCell::ChangeVolumeExponentialRelaxationEquation(
    real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus,
    real_t relaxation_rate_fluid) {
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  // Exponential relaxation towards the target volume
  const real_t current_total_volume = GetVolume();
  const real_t fluid_fraction = GetFluidFraction();
  const real_t nuclear_volume = GetNuclearVolume();

  const real_t current_nuclear_solid = nuclear_volume * (1 - fluid_fraction);
  const real_t current_cytoplasm_solid =
      (current_total_volume - nuclear_volume) * (1 - fluid_fraction);

  const real_t current_fluid = fluid_fraction * current_total_volume;

  // Update fluid volume
  real_t new_fluid =
      current_fluid +
      sparams->dt_cycle * relaxation_rate_fluid *
          (GetTargetFractionFluid() * current_total_volume - current_fluid);
  // Clamp to zero to prevent negative volumes
  if (new_fluid < 0.0) {
    new_fluid = 0.0;
  }

  const real_t nuclear_fluid =
      new_fluid * (nuclear_volume / current_total_volume);
  // real_t cytoplasm_fluid = new_fluid - nuclear_fluid;

  real_t nuclear_solid = current_nuclear_solid +
                         sparams->dt_cycle * relaxation_rate_nucleus *
                             (GetTargetNucleusSolid() - current_nuclear_solid);
  // Clamp to zero to prevent negative volumes
  if (nuclear_solid < 0.0) {
    nuclear_solid = 0.0;
  }

  const real_t target_cytoplasm_solid =
      GetTargetRelationCytoplasmNucleus() * GetTargetNucleusSolid();
  real_t cytoplasm_solid =
      current_cytoplasm_solid +
      sparams->dt_cycle * relaxation_rate_cytoplasm *
          (target_cytoplasm_solid - current_cytoplasm_solid);
  // Clamp to zero to prevent negative volumes
  if (cytoplasm_solid < 0.0) {
    cytoplasm_solid = 0.0;
  }

  const real_t new_total_solid = nuclear_solid + cytoplasm_solid;

  const real_t total_nuclear = nuclear_solid + nuclear_fluid;

  // real_t total_cytoplasm= cytoplasm_solid + cytoplasm_fluid;

  const real_t new_volume = new_total_solid + new_fluid;

  // Avoid division by zero
  const real_t new_fraction_fluid = new_fluid / (kEpsilon + new_volume);

  // Update the cell's properties
  // if the volume has changed
  if (new_volume != current_total_volume) {
    SetVolume(new_volume);
    // Update constants for all ConsumptionSecretion of Oxygen, Glucose and
    // Immunostimulatory Factors
    ComputeConstantsConsumptionSecretion();
  }
  SetFluidFraction(new_fraction_fluid);
  SetNuclearVolume(total_nuclear);
}

// compute Displacement
Real3 TumorCell::CalculateDisplacement(const InteractionForce* force,
                                       real_t squared_radius, real_t /*dt*/) {
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a
  // cleaner code instead of hardcoding it here
  squared_radius = sparams->squared_max_distance_neighbors_force;

  Real3 translation_velocity_on_point_mass{0, 0, 0};

  // We check for every neighbor object if they touch us, i.e. push us
  // away and agreagate the velocities

  uint64_t non_zero_neighbor_forces = 0;
  if (!IsStatic()) {
    ExecutionContext* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto calculate_neighbor_forces =
        L2F([&](Agent* neighbor, real_t /*squared_distance*/) {
          Real4 neighbor_force = force->Calculate(this, neighbor);
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
  // position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 *
  // velocity(t - dt) ]
  movement_at_next_step += translation_velocity_on_point_mass * sparams->dnew +
                           older_velocity_ * sparams->dold;

  older_velocity_ = translation_velocity_on_point_mass;

  // Clamp the movement if it surpasses the z boundaries.
  Real3 current_position = GetPosition();
  const double current_z = current_position[2];
  double& movement_z = movement_at_next_step[2];
  const double min_z = sparams->bounded_space_min_allowed_z;
  const double max_z = sparams->bounded_space_max_allowed_z;
  const double next_z = current_z + movement_z;
  if (next_z < min_z) {
      movement_z = min_z - current_z;
  } else if (next_z > max_z) {
      movement_z = max_z - current_z;
  }

  // Displacement
  return movement_at_next_step;
}

// Compute new oxygen or immunostimulatory factor concentration after
// consumption/ secretion
real_t TumorCell::ConsumeSecreteSubstance(int substance_id,
                                          real_t old_concentration) {
  real_t res = 0.0;
  if (substance_id == oxygen_dgrid_->GetContinuumId()) {
    // consuming oxygen
    res = (old_concentration + constant1_oxygen_) / constant2_oxygen_;
  } else if (IsImmunostimulatoryFactorDefined() && substance_id == immunostimulatory_factor_dgrid_->GetContinuumId()) {
    // There is a definde immunostimulatory factor grid and this is the one being updated
    // secreting immunostimulatory factor
    res = (old_concentration + constant1_immunostimulatory_factor_) /
          constant2_immunostimulatory_factor_;
  } else if (IsGlucoseDefined() && substance_id == glucose_dgrid_->GetContinuumId()) {
    // There is a defined glucose grid and this is the one being updated
    // CAR-T do not change glucose levels
    res = (old_concentration + constant1_glucose_) / constant2_glucose_;
  } else {
    throw std::invalid_argument("Unknown substance id: " +
                                std::to_string(substance_id));
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
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  const real_t new_volume = GetVolume();
  // compute the constants for the differential equation explicit solution: for
  // oxygen glucose and immunostimulatory factor
  // dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt
  // · (V_k / V_voxel) · S_k · ρ*_k)
  constant1_oxygen_ = 0.;
  // 1 + dt*(cell_volume/voxel_volume)*(quantity_secretion +
  // quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  //  Scale by the volume of the cell in the Voxel and time step
  constant2_oxygen_ = 1 + sparams->dt_substances *
                              (new_volume / sparams->voxel_volume) *
                              (oxygen_consumption_rate_);
  if (IsImmunostimulatoryFactorDefined()) {
    // dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt
    // · (V_k / V_voxel) · S_k · ρ*_k)
    // Scale by the volume of the cell in the Voxel and time step
    constant1_immunostimulatory_factor_ =
        immunostimulatory_factor_secretion_rate_ *
        sparams->saturation_density_immunostimulatory_factor *
        sparams->dt_substances * (new_volume / sparams->voxel_volume);
    // 1 + dt*(cell_volume/voxel_volume)*(quantity_secretion +
    // quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
    // Scale by the volume of the cell in the Voxel and time step
    constant2_immunostimulatory_factor_ =
        1 + sparams->dt_substances * (new_volume / sparams->voxel_volume) *
                (immunostimulatory_factor_secretion_rate_);
  }

  if (IsGlucoseDefined()) {
    // dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt
    // · (V_k / V_voxel) · S_k · ρ*_k)
    constant1_glucose_ = 0.;
    // 1 + dt*(cell_volume/voxel_volume)*(quantity_secretion +
    // quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
    constant2_glucose_ = 1 + sparams->dt_substances *
                                (new_volume / sparams->voxel_volume) *
                                (glucose_consumption_rate_);
  }
}

void TumorCell::StartApoptosis() {
  // If the cell is already dead, do nothing
  if (IsDead()) {
    return;
  }
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();

  // The cell Dies
  SetState(TumorCellState::kApoptotic);

  // Reset timer_state
  SetTimerState(0);
  // Set type to 5 to indicate dead cell
  SetType(TumorCellType::kType5);
  // Set target volume to 0 (the cell shrinks)
  SetTargetCytoplasmSolid(0.0);
  SetTargetNucleusSolid(0.0);
  SetTargetFractionFluid(0.0);
  SetTargetRelationCytoplasmNucleus(0.0);
  // Reduce oxygen consumption
  SetOxygenConsumptionRate(GetOxygenConsumptionRate() *
                           sparams->reduction_consumption_dead_cells);
  // Stop Immunostimulatory Factor Secretion
  SetImmunostimulatoryFactorSecretionRate(0.0);
  // Update constants for consumption/secretion differential equation solving
  ComputeConstantsConsumptionSecretion();
}

/// Main behavior executed at each simulation step
void StateControlGrowProliferate::Run(Agent* agent) {
  Simulation* sim = Simulation::GetActive();
  const auto* sparams = sim->GetParam()->Get<SimParam>();
  // Run only every sparams->dt_cycle minutes, fmod does not work with the type
  // returned by GetSimulatedTime()
  if (sim->GetScheduler()->GetSimulatedSteps() %
          sparams->steps_per_cell_cycle !=
      0) {
    return;
  }

  if (auto* cell = dynamic_cast<TumorCell*>(agent)) {
    if (cell->IsAttachedToCart()) {
      // If the cell is attached to a cart, skip the state control and growth
      return;
    }
    switch (cell->GetState()) {
      case TumorCellState::kAlive: {
        // the cell is growing to real_t its size before mitosis
        // Increase timer_state to track time in this state (sparams->dt_cycle
        // minutes per step)
        cell->SetTimerState(cell->GetTimerState() + sparams->dt_cycle);

          // Oxygen levels
          const Real3 current_position = cell->GetPosition();
          // Pointer to the oxygen diffusion grid
          DiffusionGrid* oxygen_dgrid = cell->GetOxygenDiffusionGrid();
          const real_t oxygen_level = oxygen_dgrid->GetValue(current_position);

          // Glucose levels if a glucose diffusion grid is defined
          real_t glucose_level = 0.0;
          if (cell->IsGlucoseDefined()) {
            DiffusionGrid* glucose_dgrid = cell->GetGlucoseDiffusionGrid();
            glucose_level = glucose_dgrid->GetValue(current_position);
          }

        // Enter necrosis if oxygen or glucose levels are too low
        if (ShouldEnterNecrosis(oxygen_level, glucose_level, cell)) {
          // Exit the function to prevent further processing
          return;
        }
        ManageLivingCell(cell, oxygen_level, glucose_level);
        break;
      }
      case TumorCellState::kNecroticSwelling: {
        // the cell is swelling before lysing
        // Increase timer_state to track time in this state (sparams->dt_cycle
        // minutes per step)
        cell->SetTimerState(cell->GetTimerState() + sparams->dt_cycle);
        // volume change
        //  The cell swells
        cell->ChangeVolumeExponentialRelaxationEquation(
            sparams
                ->volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell,
            sparams
                ->volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell,
            sparams->volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell);

        // If the cell has swollen to 2 times its original volume, it lyses
        if (cell->GetVolume() >= 2 * sparams->default_volume_new_tumor_cell) {
          // Change state to necrotic lysed
          cell->SetState(TumorCellState::kNecroticLysed);
          // Reset timer_state
          cell->SetTimerState(0);
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
          // Update constants for all ConsumptionSecretion of Oxygen and
          // Immunostimulatory Factors
          cell->ComputeConstantsConsumptionSecretion();
        }
        break;
      }
      case TumorCellState::kNecroticLysed: {
        // the cell is shrinking and will be removed after a certain time
        // Increase timer_state to track time in this state (sparams->dt_cycle
        // minutes per step)
        cell->SetTimerState(cell->GetTimerState() + sparams->dt_cycle);
        // volume change
        //  The cell shrinks
        cell->ChangeVolumeExponentialRelaxationEquation(
            sparams->volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell,
            sparams->volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell,
            sparams->volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell);
        // If the timer_state exceeds the time to transition (this is a fixed
        // duration transition)
        if (sparams->time_lysis < cell->GetTimerState()) {
          // remove the cell from the simulation
          ExecutionContext* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      case TumorCellState::kApoptotic: {
        // Increase timer_state to track time in this state (sparams->dt_cycle
        // minutes per step)
        cell->SetTimerState(cell->GetTimerState() + sparams->dt_cycle);

        //  The cell shrinks
        cell->ChangeVolumeExponentialRelaxationEquation(
            sparams->volume_relaxation_rate_cytoplasm_apoptotic_cells,
            sparams->volume_relaxation_rate_nucleus_apoptotic_cells,
            sparams->volume_relaxation_rate_fluid_apoptotic_cells);
        // If the timer_state exceeds the time to transition (this is a fixed
        // duration transition)
        if (sparams->time_apoptosis < cell->GetTimerState()) {
          // remove the cell from the simulation
          ExecutionContext* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      default: {
        Log::Error("StateControlGrowProliferate::Run",
                   "Unknown TumorCellState");
        break;
      }
    }
  } else {
    Log::Error("StateControlGrowProliferate::Run",
               "SimObject is not a TumorCell");
  }
}

// ManageLivingCell function to handle living cell behavior
void StateControlGrowProliferate::ManageLivingCell(TumorCell* cell,
                                                   real_t oxygen_level, real_t glucose_level) {

  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  // volume change
  // If this tumoral cell is already bigger than the target volume, the volume does not change since it is already big enough to divide and it should not shrink
  if(cell->GetVolume() < cell->GetTargetTotalVolume()) {
    // If there is a glucose diffusion grid defined it affects the speed of growth of the cell
    // Initialize multiplier
    real_t multiplier_glucose = 1.0;
    if(cell->IsGlucoseDefined()) {
      //if there is a definde glucose gradient in the simulation
      // glucose threshold for considering an effect on the growth speed
      if (glucose_level < sparams->glucose_saturation_for_tumor_cell_growth) {
        multiplier_glucose = (glucose_level - sparams->glucose_limit_for_tumor_cell_growth) /
                    (sparams->glucose_saturation_for_tumor_cell_growth -
                      sparams->glucose_limit_for_tumor_cell_growth);
      }
      // If glucose is below the limit, set multiplier to 0
      if (glucose_level < sparams->glucose_limit_for_tumor_cell_growth) {
        multiplier_glucose = 0.0;
      }
    }
    //change volume change speed depending on the vailable glucose. This also affects ploriferation because if the cells do not grow big enough they cannot divide
    cell->ChangeVolumeExponentialRelaxationEquation(
        sparams->volume_relaxation_rate_alive_tumor_cell_cytoplasm * multiplier_glucose,
        sparams->volume_relaxation_rate_alive_tumor_cell_nucleus * multiplier_glucose,
        sparams->volume_relaxation_rate_alive_tumor_cell_fluid * multiplier_glucose);
  }

  // cell state control
  // The division rate depends on the available resources
  // Oxygen
  // Initialize multiplier
  real_t multiplier_oxygen = 1.0;
  // oxygen threshold for considering an effect on the proliferation cycle
  if (oxygen_level < sparams->oxygen_saturation_for_proliferation) {
    multiplier_oxygen = (oxygen_level - sparams->oxygen_limit_for_proliferation) /
                 (sparams->oxygen_saturation_for_proliferation -
                  sparams->oxygen_limit_for_proliferation);
  }
  // If oxygen is below the limit, set multiplier to 0
  if (oxygen_level < sparams->oxygen_limit_for_proliferation) {
    multiplier_oxygen = 0.0;
  }
  // Calculate the rate of state change based on oxygen level and oncoprotein
  // (min^-1)
  real_t final_rate_transition = cell->GetTransformationRandomRate() *
                                       multiplier_oxygen * cell->GetOncoproteinLevel();

  // Check the dimensions of the cell, if it is too small, it cannot divide: it needs to grow more before dividing
  if(cell->GetVolume() < cell->GetTargetTotalVolume()*sparams->minimum_tumor_cell_target_volume_fraction_for_division) {
    final_rate_transition=0;
  }

  // Calculate the time to wait (in minutes)
  real_t time_to_wait = kTimeTooLarge;
  if (final_rate_transition > 0) {
    // Calculate the time to transition (in minutes)
    time_to_wait = 1. / final_rate_transition;
  }
  // If the timer_state exceeds the time to transition (this is a fixed duration
  // transition)
  if (time_to_wait < cell->GetTimerState()) {
    // mitosis: cell divides
    cell->SetState(TumorCellState::kAlive);
    cell->Divide();
    // Reset timer_state
    cell->SetTimerState(0);
  }
}

// computes the probability of the cell entering necrosis
bool StateControlGrowProliferate::ShouldEnterNecrosis(real_t oxygen_level, real_t glucose_level,
                                                      TumorCell* cell) {
  Simulation* sim = Simulation::GetActive();
  const auto* sparams = sim->GetParam()->Get<SimParam>();
  // necrosis probability
  // Oxygen
  // Default multiplier for necrosis probability because of low oxygen level
  real_t multiplier = 0.0;
  // oxygen threshold for considering necrosis
  if (oxygen_level < sparams->oxygen_limit_for_necrosis) {
    multiplier = (sparams->oxygen_limit_for_necrosis - oxygen_level) /
                 (sparams->oxygen_limit_for_necrosis -
                  sparams->oxygen_limit_for_necrosis_maximum);
  }
  // threshold for maximum necrosis probability
  if (oxygen_level < sparams->oxygen_limit_for_necrosis_maximum) {
    multiplier = 1.0;
  }
  // Calculate the probability of necrosis based on oxygen level
  const real_t maximum_necrosis_rate_oxygen_multiplier = sparams->maximum_necrosis_lack_of_oxygen_rate * multiplier;
  // Glucose (only if glucose is defined in the simulation)
  // Default multiplier for necrosis probability because of low glucose level
  multiplier = 0.0;
  if (cell->IsGlucoseDefined()) {
    //there is a glucose diffusion grid defined in the simulation
    // glucose threshold for considering necrosis
    if (glucose_level < sparams->glucose_limit_for_necrosis) {
      multiplier = (sparams->glucose_limit_for_necrosis - glucose_level) /
                  (sparams->glucose_limit_for_necrosis -
                    sparams->glucose_limit_for_necrosis_maximum);
    }
    // threshold for maximum necrosis probability
    if (glucose_level < sparams->glucose_limit_for_necrosis_maximum) {
      multiplier = 1.0;
    }
  }
  // Calculate the probability of necrosis based on glucose level. If no glucose gradient is defined it will be zero
  const real_t maximum_necrosis_rate_glucose_multiplier = sparams->maximum_necrosis_lack_of_glucose_rate * multiplier;
  // Random natural causes
  const real_t current_basal_death_probability = cell->GetBasalDeathProbability();
  // Final probability: multiply by sparams->dt_cycle since each timestep is sparams->dt_cycle minutes
  // The probability of necrosis is calculated as the complement of the product of the complements of the individual probabilities, scaled by the time step.
  const real_t probability_necrosis =
      sparams->dt_cycle * (1.0 - (1.0 - maximum_necrosis_rate_oxygen_multiplier) *
                                  (1.0 - maximum_necrosis_rate_glucose_multiplier) *
                                  (1.0 - current_basal_death_probability));

  Random* random = sim->GetRandom();
  const bool enter_necrosis = random->Uniform(0, 1) < probability_necrosis;
  // If the random number is less than the probability, enter necrosis
  if (enter_necrosis) {
    // If oxygen is too low, enter necrosis
    cell->SetState(TumorCellState::kNecroticSwelling);
    // Reset timer_state
    cell->SetTimerState(0);

    // Stop Secretion and reduce consumption
    //  Stop secretion
    cell->SetImmunostimulatoryFactorSecretionRate(0.0);
    // Reduce consumption
    cell->SetOxygenConsumptionRate(cell->GetOxygenConsumptionRate() *
                                   sparams->reduction_consumption_dead_cells);
    cell->SetGlucoseConsumptionRate(cell->GetGlucoseConsumptionRate() *
                                   sparams->reduction_consumption_dead_cells);
    // Update constants for all ConsumptionSecretion of Oxygen and
    // Immunostimulatory Factors
    cell->ComputeConstantsConsumptionSecretion();

    // The cell will swell getting filled with fluid
    cell->SetTargetCytoplasmSolid(0);
    cell->SetTargetNucleusSolid(0);
    // Set target fraction of fluid to 1.0
    cell->SetTargetFractionFluid(1.0);
    cell->SetTargetRelationCytoplasmNucleus(0.0);
    // Set type to 5 to indicate dead cell
    cell->SetType(TumorCellType::kType5);
  }
  return enter_necrosis;  // Return whether the cell entered necrosis
}

}  // namespace bdm
