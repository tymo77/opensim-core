#ifndef MOCO_MOCOENERGYGOAL_H
#define MOCO_MOCOENERGYGOAL_H
/* -------------------------------------------------------------------------- *
 * OpenSim Moco: MocoEnergyGoal.h                                             *
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

#include "MocoGoal.h"

namespace OpenSim {

/// Create pair of variables for use with a MocoPowerConstraint.
class OSIMMOCO_API MocoEnergyGoalPair : public Object {
    OpenSim_DECLARE_CONCRETE_OBJECT(MocoEnergyGoalPair, Object);

public:
    OpenSim_DECLARE_PROPERTY(
        torque_var,
        std::string,
        "The torque control variable for the pair.");
    OpenSim_DECLARE_PROPERTY(
        speed_var,
        std::string,
        "The speed state variable for the pair.");
    OpenSim_DECLARE_PROPERTY(
        negate,
        bool,
        "Should the first and second variables have opposite signs? Defaults to false.");


    MocoEnergyGoalPair();
    MocoEnergyGoalPair(std::string torque_var, std::string speed_var, bool negate = false);

private:
    void constructProperties();
};

class OSIMMOCO_API MocoEnergyGoal : public MocoGoal {
    OpenSim_DECLARE_CONCRETE_OBJECT(
        MocoEnergyGoal, MocoGoal);

public:
    MocoEnergyGoal();
    MocoEnergyGoal(std::string name) : MocoGoal(std::move(name)) {
        constructProperties();
    }
    MocoEnergyGoal(std::string name, double weight)
        : MocoGoal(std::move(name), weight) {
        constructProperties();
    }

    OpenSim_DECLARE_PROPERTY(
        smooth_scale,
        double,
        "The inverse order of magnitude for the region where smoothing occurs, in power units. Default: 1.0.");
    OpenSim_DECLARE_PROPERTY(
        regen_efficiency,
        double,
        "Efficiency of regeneration. Should be 0 to 1. Default: 0.0."
    );

    void addPair(MocoEnergyGoalPair pair) {
        append_pairs(std::move(pair));
    }
    void addNegatedPair(MocoEnergyGoalPair pair) {
        pair.set_negate(true);
        append_pairs(std::move(pair));
    }

    void setRegenEfficiency(double r) { set_regen_efficiency(r); }
    double getRegenEfficiency() const { return get_regen_efficiency(); }

    void setSmoothScale(double s) { set_smooth_scale(s); }
    double getSmoothScale() const { return get_smooth_scale(); }

    /// Set if the goal should be divided by the displacement of the system's
        /// center of mass over the phase.
    void setDivideByDisplacement(bool tf) { set_divide_by_displacement(tf); }
    bool getDivideByDisplacement() const {
        return get_divide_by_displacement();
    }

    /** Set if the goal should be divided by the duration the phase. */
    void setDivideByDuration(bool tf) { set_divide_by_duration(tf); }
    bool getDivideByDuration() const { return get_divide_by_duration(); }

protected:
    Mode getDefaultModeImpl() const override { return Mode::Cost; }
    bool getSupportsEndpointConstraintImpl() const override { return true; }
    void initializeOnModelImpl(const Model&) const override;
    void calcIntegrandImpl(
            const IntegrandInput& input, SimTK::Real& integrand) const override;
    void calcGoalImpl(
            const GoalInput& input, SimTK::Vector& cost) const override;
    void printDescriptionImpl() const override;

private:
    OpenSim_DECLARE_LIST_PROPERTY(pairs, MocoEnergyGoalPair,
        "Pairs of speeds and torques for calculating and limiting mechanical power.");
    OpenSim_DECLARE_PROPERTY(divide_by_displacement, bool,
        "Divide by the model's displacement over the phase (default: "
            "false)");
    OpenSim_DECLARE_PROPERTY(
            divide_by_duration, bool, "Divide by the total time of the phase.");
    void constructProperties();
    mutable std::vector<std::tuple<int, int, int>> m_indices;
    mutable std::vector<std::pair<std::string, std::string>> m_names;
    static double softplus(const double x, const double k) {
        double kx = k * x;
        if (kx > 700.0) {
            // If the multiple kx is too large, we can get overflow.
            // So here is a really messy rule past that point.
            return x;
        }
        else {
            return log(1 + exp(kx)) / k;
        }
    };
};

} // namespace OpenSim

#endif // MOCO_MOCOENERGYGOAL_H
