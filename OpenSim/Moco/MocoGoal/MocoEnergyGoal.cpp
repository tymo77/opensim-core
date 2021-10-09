/* -------------------------------------------------------------------------- *
 * OpenSim Moco: MocoEnergyGoal.cpp                                           *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2019 Stanford University and the Authors                     *
 *                                                                            *
 * Author(s): Christopher Dembia                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0          *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "MocoEnergyGoal.h"

using namespace OpenSim;

//=============================================================================
//  MocoEnergyGoalPair
//=============================================================================

MocoEnergyGoalPair::MocoEnergyGoalPair() { constructProperties(); }

MocoEnergyGoalPair::MocoEnergyGoalPair(
    std::string torque_var, std::string speed_var, bool negate) {
    constructProperties();
    set_torque_var(torque_var);
    set_speed_var(speed_var);
    set_negate(negate);
}

void MocoEnergyGoalPair::constructProperties() {
    constructProperty_torque_var("");
    constructProperty_speed_var("");
    constructProperty_negate(false);
}

//=============================================================================
//  MocoEnergyGoal
//=============================================================================

MocoEnergyGoal::MocoEnergyGoal() {
    constructProperties();
}

void MocoEnergyGoal::constructProperties() {
    constructProperty_pairs();
    constructProperty_divide_by_displacement(false);
    constructProperty_divide_by_duration(false);
    constructProperty_smooth_scale(1.0);
    constructProperty_regen_efficiency(0.0);
}

void MocoEnergyGoal::initializeOnModelImpl(
        const Model& model) const {

    // Get all expected control names.
    auto controlNames = createControlNamesFromModel(model);

    // Check that the model controls are in the correct order.
    checkOrderSystemControls(model);

    // Get the state and control index maps.
    auto allSysYIndices =        createSystemYIndexMap(model);
    auto systemControlIndexMap = createSystemControlIndexMap(model);

    int nPairs = getProperty_pairs().size();
    for (int i = 0; i < nPairs; ++i) {
        const auto& torque_path = get_pairs(i).get_torque_var();
        if (systemControlIndexMap.count(torque_path) == 0) {
            OPENSIM_THROW_FRMOBJ(
                    Exception, "Unrecognized control '" + torque_path + "'.");
        }

        const auto& speed_path = get_pairs(i).get_speed_var();
        if (allSysYIndices.count(speed_path) == 0) {
            OPENSIM_THROW_FRMOBJ(
                    Exception, "Unrecognized state '" + speed_path + "'.");
        }

        int controlIndex = systemControlIndexMap[torque_path];
        int stateIndex = allSysYIndices[speed_path];
        m_names.emplace_back(torque_path, speed_path);
        m_indices.emplace_back(controlIndex, stateIndex,
            get_pairs(i).get_negate() ? -1 : 1);
    }

    /*setRequirements(1, 1,
            get_divide_by_displacement() ? SimTK::Stage::Position
                                         : SimTK::Stage::Model);*/
    setRequirements(1, 1, SimTK::Stage::Time);
        
}

void MocoEnergyGoal::calcIntegrandImpl(
        const IntegrandInput& input, SimTK::Real& integrand) const {
    
    auto& state = input.state;

    // Iterate through the pairs.
    const auto& s = state.getY();
    const auto& u = input.controls;
    const auto& t = input.time;

    double ss = get_smooth_scale();
    double eff = get_regen_efficiency();

    int idx_control, idx_state;
    double torque, speed, sign, power;
    for (const auto& indices : m_indices) {

        idx_control = std::get<0>(indices);
        idx_state   = std::get<1>(indices);

        torque = u[idx_control];
        speed  = s[idx_state];
        sign = std::get<2>(indices);

        power = torque * speed * sign;

        // The positive power is softplused to smoothly transition from slope of 1 to zero.
        double positive = softplus(power, ss);

        // The negative power is softplused to smoothly transition from slope of eff (regen efficiency) to zero.
        double negative = -eff * softplus(-power, ss);

        integrand += (positive + negative);
    }
}

void MocoEnergyGoal::calcGoalImpl(
        const GoalInput& input, SimTK::Vector& cost) const {
    cost[0] = input.integral;
    if (get_divide_by_displacement()) {
        cost[0] /=
                calcSystemDisplacement(input.initial_state, input.final_state);
    }

    if (get_divide_by_duration()) {
        cost[0] /= (input.final_time - input.initial_time);
    }
}

void MocoEnergyGoal::printDescriptionImpl() const {
    int nPairs = getProperty_pairs().size();
    for (int i = 0; i < nPairs; ++i) {
        const auto& torque_path = get_pairs(i).get_torque_var();
        const auto& speed_path = get_pairs(i).get_speed_var();
        const auto& negate = get_pairs(i).get_negate();
        log_cout("        control: {}, state: {}, negate: {}",
            torque_path, speed_path, negate);
    }
}