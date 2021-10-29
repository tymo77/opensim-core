/* -------------------------------------------------------------------------- *
 * OpenSim Moco: example2DWalking.cpp                                         *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2017-19 Stanford University and the Authors                  *
 *                                                                            *
 * Author(s): Antoine Falisse                                                 *
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

/// This example features two different optimal control problems:
///  - The first problem is a tracking simulation of walking.
///  - The second problem is a predictive simulation of walking.
///
/// The code is inspired from Falisse A, Serrancoli G, Dembia C, Gillis J,
/// De Groote F: Algorithmic differentiation improves the computational
/// efficiency of OpenSim-based trajectory optimization of human movement.
/// PLOS One, 2019.
///
/// Model
/// -----
/// The model described in the file '2D_gait.osim' included in this file is a
/// modified version of the 'gait10dof18musc.osim' available within OpenSim. We
/// replaced the moving knee flexion axis by a fixed flexion axis, replaced the
/// Millard2012EquilibriumMuscles by DeGrooteFregly2016Muscles, and added
/// SmoothSphereHalfSpaceForces (two contacts per foot) to model the
/// contact interactions between the feet and the ground.
///
/// Do not use this model for research. The path of the gastroc muscle contains
/// an error--the path does not cross the knee joint.
///
/// Data
/// ----
/// The coordinate data included in the 'referenceCoordinates.sto' comes from
/// predictive simulations generated in Falisse et al. 2019.

#include <OpenSim/Actuators/DCMotor.h>
#include <OpenSim/Actuators/TorqueActuator.h>
#include <OpenSim/Actuators/SpringGeneralizedForce.h>
#include <OpenSim/Common/STOFileAdapter.h>
#include <OpenSim/Common/TableUtilities.h>
#include <OpenSim/Moco/osimMoco.h>
#include <cmath>

using namespace OpenSim;

// Helper function for computing the intersection of two string vectors.
std::vector<std::string> setIntersection(
        std::vector<std::string> v1, std::vector<std::string> v2) {
    std::vector<std::string> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(
            v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
    return v3;
}

std::vector<double> boundIntersection(
        std::vector<double> b1, std::vector<double> b2) {

    // Ensure the two bound sets are sorted.
    std::sort(b1.begin(), b1.end());
    std::sort(b2.begin(), b2.end());

    // Compare the lower bounds.
    if (b1[0] < b2[0]) { b1[0] = b2[0]; } // Lower bound of b2 is greater than lower bound of b1.
    if (b1[1] > b2[1]) { b1[1] = b2[1]; } // Upper bound of b2 is less than the upper bound of b2.
    return b1;
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// //
// ===============================================================================
// //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
MocoTrajectory genInitialGuess(MocoStudy study, MocoTrajectory initial_guess) {

    // Fix the problem of somtimes NaN's in slacks by just copying the last
    // value. I'm not sure why built-in reshape and interpolate aren't catching
    // this.
    // TODO: Interpolate states and controls with appropriate degree.
    auto& slacks = initial_guess.getSlacksTrajectory();
    int N_rows = initial_guess.getNumTimes();
    int N_cols = initial_guess.getSlackNames().size();
    for (int i_col = 0; i_col < N_cols; i_col++) {
        for (int i_row = 1; i_row < N_rows; i_row++) {
            double val = slacks[i_row][i_col];
            if (SimTK::isNaN(val)) {
                slacks[i_row][i_col] = slacks[i_row - 1][i_col];
            }
        }
    }

    // Use the coordinate trajectory initial guess as best as possible to form
    // the guess.
    MocoCasADiSolver& solver = study.initCasADiSolver();
    MocoTrajectory generated_guess = solver.createGuess("bounds");

    // Guess copied and interpolated from the provided initial guess.
    int Nstates_guess = generated_guess.getNumTimes();

    // Resample the provided guess to be the same length as the generated
    // guess.
    SimTK::Vector t_vec = SimTK::Vector(Nstates_guess);
    double t0 = initial_guess.getInitialTime();
    double tf = initial_guess.getFinalTime();
    double dt = (tf - t0) / ((double)Nstates_guess - 1);
    for (int i = 0; i < Nstates_guess; i++) { t_vec[i] = t0 + dt * i; }

    t_vec[0] = initial_guess.getInitialTime();
    t_vec[Nstates_guess - 1] = initial_guess.getFinalTime();
    initial_guess.resample(t_vec);

    const auto& common_states = setIntersection(
            initial_guess.getStateNames(), generated_guess.getStateNames());
    const auto& common_controls = setIntersection(
            initial_guess.getControlNames(), generated_guess.getControlNames());
    const auto& common_derivatives =
            setIntersection(initial_guess.getDerivativeNames(),
                    generated_guess.getDerivativeNames());
    const auto& common_multipliers =
            setIntersection(initial_guess.getMultiplierNames(),
                    generated_guess.getMultiplierNames());
    const auto& common_slacks = setIntersection(
            initial_guess.getSlackNames(), generated_guess.getSlackNames());

    // Set same time.
    generated_guess.setTime(initial_guess.getTime());

    // Copy all the common states from the provided initial guess if
    // possible.
    std::cout << "The following states were found in the example:" << std::endl;
    for (const auto& s : common_states) {
        std::cout << s << std::endl;
        generated_guess.setState(s, initial_guess.getState(s));
    }

    // Copy all the common controls from the provided initial guess if
    // possible.
    std::cout << "The following controls were found in the example:"
              << std::endl;
    for (const auto& s : common_controls) {
        std::cout << s << std::endl;
        generated_guess.setControl(s, initial_guess.getControl(s));
    }

    // Copy all the common derivatives from the provided initial guess if
    // possible.
    std::cout << "The following derivatives were found in the example:"
              << std::endl;
    for (const auto& s : common_derivatives) {
        std::cout << s << std::endl;
        generated_guess.setDerivative(s, initial_guess.getDerivative(s));
    }

    // Copy all the common multipliers from the provided initial guess if
    // possible.
    std::cout << "The following multipliers were found in the example:"
              << std::endl;
    for (const auto& s : common_multipliers) {
        std::cout << s << std::endl;
        generated_guess.setMultiplier(s, initial_guess.getMultiplier(s));
    }

    // Copy all the common slacks from the provided initial guess if
    // possible.
    std::cout << "The following slacks were found in the example:" << std::endl;
    for (const auto& s : common_slacks) {
        std::cout << s << std::endl;
        generated_guess.setSlack(s, initial_guess.getSlack(s));
    }

    //generated_guess.write(guess_fn + "_guess.sto");
    return generated_guess;
}

std::vector<double> linkageKinematicsT2(
        double t2, double l0, double l1, double l2, double l3) {
    using SimTK::Pi;

    // Solving for theta 3.
    double l12sq = pow(l1, 2) + pow(l2, 2) - 2 * l1 * l2 * cos(Pi - t2);
    double cosbeta = (l12sq - pow(l0, 2) - pow(l3, 2)) / (-2 * l0 * l3);
    double sinbeta = sqrt(1.0 - pow(cosbeta, 2));
    double beta = atan2(sinbeta, cosbeta);
    double t3 = Pi / 2 - beta;

    // Solving for theta 1.
    double c3 = cos(t3);
    double s3 = sin(t3);
    double c2 = cos(t2);
    double s2 = sin(t2);
    double D = pow(l1, 2) + pow(l2, 2) + 2 * l1 * l2 * c2;

    double c1 = (l3 * c3 * (l1 + l2 * c2) + l2 * s2 * (l3 * s3 - l0)) / D;
    double s1 = (-l2 * s2 * l3 * c3 + (l1 + l2 * c2) * (l3 * s3 - l0)) / D;
    double t1 = atan2(s1, c1);

    return {t1, t2, t3};
}

std::vector<std::vector<double>> linkageKinematicsT3(
        double t3, double l0, double l1, double l2, double l3) {
    using SimTK::Pi;
    
    double s3 = sin(t3);
    double c3 = cos(t3);

    double A = -2 * l1 * l3 * c3;
    double B = 2 * l0 * l1 - 2 * l1 * l3 * s3;
    double C = pow(l0, 2) + pow(l1, 2) + pow(l3, 2) - pow(l2, 2) -
               2 * l0 * l3 * s3;
    double D = pow(B, 2) - pow(C, 2) + pow(A, 2);
    double tsol1, tsol2;
    if (D < 0) {
        // In this case, we find the "reaching angle".
        tsol1 = (-B ) / (C - A);
        tsol2 = (-B ) / (C - A);
    } else {
        tsol1 = (-B + sqrt(D)) / (C - A);
        tsol2 = (-B - sqrt(D)) / (C - A);
    }

    double t1_sol1 = 2 * atan(tsol1);
    double t1_sol2 = 2 * atan(tsol2);

    double c12_sol1 = (l3 * c3 - l1 * cos(t1_sol1)) / l2;
    double s12_sol1 = (l3 * s3 - l0 - l1 * sin(t1_sol1)) / l2;
    double c12_sol2 = (l3 * c3 - l1 * cos(t1_sol2)) / l2;
    double s12_sol2 = (l3 * s3 - l0 - l1 * sin(t1_sol2)) / l2;

    double t2_sol1 = atan2(s12_sol1, c12_sol1) - t1_sol1;
    double t2_sol2 = atan2(s12_sol2, c12_sol2) - t1_sol2;

    return {{t1_sol1, t2_sol1, t3}, {t1_sol2, t2_sol2, t3}};
}

void configureLinkage(double l0, double l1, double l2, double l3, Model& model, double ub, double lb) {

    // Get the angles for t1 and t3 that correspont to the upper and lower bounds of t2.
    auto angles1 = linkageKinematicsT2(lb, l0, l1, l2, l3);
    auto angles2 = linkageKinematicsT2(ub, l0, l1, l2, l3);

    std::vector<double> t1_kinbounds = {angles1[0], angles2[0]};
    std::vector<double> t2_kinbounds = {lb, ub}; 
    std::vector<double> t3_kinbounds = {angles1[2], angles2[2]};

    // Sort the other bounds which might not be ordered.
    std::sort(t1_kinbounds.begin(), t1_kinbounds.end());
    std::sort(t3_kinbounds.begin(), t3_kinbounds.end());

    std::cout << "Implied kinematic limits, t1: " << t1_kinbounds[0] << " to " << t1_kinbounds[1]
              << std::endl;
    std::cout << "Implied kinematic limits, t2: " << t2_kinbounds[0] << " to "
              << t2_kinbounds[1] << std::endl;
    std::cout << "Implied kinematic limits, t3: " << t3_kinbounds[0] << " to "
              << t3_kinbounds[1] << std::endl;

    // Now we need to go to the model and set the linkage lengths.
    model.initSystem();
    auto &s = model.initializeState();

    auto &B = model.updBodySet();
    auto &J = model.updJointSet();

    for (std::string side : {"r", "l"}) {
        auto& tibia = B.get("tibia_" + side);
        auto& talus = B.get("talus_" + side);
        auto& link1 = model.updComponent<Body>("link1_" + side);
        auto& link2 = model.updComponent<Body>("link2_" + side);

        auto& knee = J.get("knee_" + side);
        auto& ankle = J.get("ankle_" + side);
        auto& j_in = model.updComponent<Joint>("j_in_" + side);
        auto& j_passive1 = model.updComponent<Joint>("j_passive1_" + side);

        auto& pin_con =
                model.updComponent<PointOnLineConstraint>("pin_con_" + side);

        // Position of the Ankle in the Tibia frame.
        PhysicalOffsetFrame& tibia_offset_frame = knee.upd_frames(1);
        auto P_tibia_ankle =
                tibia_offset_frame.findStationLocationInAnotherFrame(
                        s, SimTK::Vec3(0, 0, 0), talus);
        double l_tibia = P_tibia_ankle[1];

        // Set the input link geometry.
        auto& link1_geo = link1.upd_attached_geometry(0);
        Brick* link1_brick = dynamic_cast<Brick*>(&(link1_geo));
        auto& link1_half_lengths = link1_brick->upd_half_lengths();
        link1_half_lengths[0] = l1 / 2;
        
        // Set the input joint transformations.
        PhysicalOffsetFrame& j_in_parent = j_in.upd_frames(0);
        PhysicalOffsetFrame& j_in_child = j_in.upd_frames(1);
        j_in_parent.set_translation({0, -l_tibia + l0, 0});
        j_in_child.set_translation({l1 / 2, 0, 0});

        // Set the coupler link geometry.
        auto& link2_geo = link2.upd_attached_geometry(0);
        Brick* link2_brick = dynamic_cast<Brick*>(&(link2_geo));
        auto& link2_half_lengths = link2_brick->upd_half_lengths();
        link2_half_lengths[0] = l2 / 2;

        // Set the passive joint transformations.
        PhysicalOffsetFrame& j_passive1_parent = j_passive1.upd_frames(0);
        PhysicalOffsetFrame& j_passive1_child = j_passive1.upd_frames(1);
        j_passive1_parent.set_translation({-l1 / 2, 0, 0});
        j_passive1_child.set_translation({l2 / 2, 0, 0});

        // Set the loop-closure constraint transformations.
        pin_con.setPointOnLine({-l2 / 2, 0, 0});
        pin_con.setPointOnFollower({l3, 0, 0});

        // Go through the model and set the coordinate bounds for each of the
        // linkage angles.
        auto& ankle_angle = model.updComponent<Coordinate>(
                "/jointset/ankle_" + side + "/ankle_angle_" + side);
        auto& motor_angle =
                model.updComponent<Coordinate>("/j_in_"  + side + "/motor_angle_" + side);
        auto& passive_angle1 = model.updComponent<Coordinate>(
                "/j_passive1_" + side + "/passive_angle1_" + side);

        // Set the default angles for the linkage angles according to the ankle
        // default.
        auto default_angles = linkageKinematicsT3(
                ankle_angle.getDefaultValue(), l0, l1, l2, l3);
        std::cout << "Default angles 1: " << default_angles[0][0] << " "
                  << default_angles[0][1] << " " << default_angles[0][2]
                  << std::endl;
        std::cout << "Default angles 1: " << default_angles[1][0] << " "
                  << default_angles[1][1] << " " << default_angles[1][2]
                  << std::endl;
        // Check which solution assembly mode is within the t2 bounds.
        int assembly_mode;
        if (default_angles[0][1] <= t2_kinbounds[1] &&
                default_angles[0][1] >= t2_kinbounds[0]) {
            assembly_mode = 0;
        } else if (default_angles[1][1] <= t2_kinbounds[1] &&
                   default_angles[1][1] >= t2_kinbounds[0]) {
            assembly_mode = 1;
        } else {
            OPENSIM_THROW(Exception, "Default configuration is not within "
                                     "perscribed bounds for passive angle");
        }
        motor_angle.setDefaultValue(default_angles[assembly_mode][0]);
        passive_angle1.setDefaultValue(default_angles[assembly_mode][1]);


        auto lower_bound_sol = linkageKinematicsT3(
                ankle_angle.getRangeMin(), l0, l1, l2, l3);
        auto upper_bound_sol = linkageKinematicsT3(
                ankle_angle.getRangeMax(), l0, l1, l2, l3);
        std::vector<double> t1_anklebounds{lower_bound_sol[assembly_mode][0],
                upper_bound_sol[assembly_mode][0]};
        std::vector<double> t2_anklebounds{lower_bound_sol[assembly_mode][1],
                upper_bound_sol[assembly_mode][1]};
        std::vector<double> t3_anklebounds{lower_bound_sol[assembly_mode][2],
                upper_bound_sol[assembly_mode][2]};
        std::sort(t1_anklebounds.begin(), t1_anklebounds.end());
        std::sort(t2_anklebounds.begin(), t2_anklebounds.end());
        std::sort(t3_anklebounds.begin(), t3_anklebounds.end());

        std::cout << "Implied limits based on " + side + " ankle motion, t1: "
                  << t1_anklebounds[0] << " to " << t1_anklebounds[1]
                  << std::endl;
        std::cout << "Implied limits based on " + side + " ankle motion, t2: "
                  << t2_anklebounds[0] << " to " << t2_anklebounds[1]
                  << std::endl;
        std::cout << "Implied limits based on " + side + " ankle motion, t3: "
                  << t3_anklebounds[0] << " to " << t3_anklebounds[1]
                  << std::endl;

        // Bounds on each joint should be the intersection of the ankle bounds
        // and the kinematic bounds.
        auto t1_finalbounds = boundIntersection(t1_kinbounds, t1_anklebounds);
        auto t2_finalbounds = boundIntersection(t2_kinbounds, t2_anklebounds);
        auto t3_finalbounds = boundIntersection(t3_kinbounds, t3_anklebounds);
        std::cout << "Final bounds, " + side + " t1: " << t1_finalbounds[0]
                  << " to "
                  << t1_finalbounds[1] << std::endl;
        std::cout << "Final bounds, " + side + " t2: " << t2_finalbounds[0]
                  << " to "
                  << t2_finalbounds[1] << std::endl;
        std::cout << "Final bounds, " + side + " t3: " << t3_finalbounds[0]
                  << " to "
                  << t3_finalbounds[1] << std::endl;

        motor_angle.setRangeMin(t1_finalbounds[0]);
        motor_angle.setRangeMax(t1_finalbounds[1]);

        passive_angle1.setRangeMin(t2_finalbounds[0]);
        passive_angle1.setRangeMax(t2_finalbounds[1]);

        ankle_angle.setRangeMin(t3_finalbounds[0]);
        ankle_angle.setRangeMax(t3_finalbounds[1]);
    }
     
}


// Fill-in the linkage angles and speeds based on the trajectory.
// Linkage angles are chosen by solving inverse kinematics.
// Speeds are determined by simply first-order finite differencing.
MocoTrajectory& generateLinkageAngles(Model& model, MocoTrajectory& mt) {
    // Symmetric coordinate values (except for pelvis_tx) and speeds.
    model.initSystem();

   auto& s = model.initializeState();

   int Ntime = mt.getNumTimes();

   // Loop through every instant in time.
   // Solve the kinematics. 
   for (int i_t = 0; i_t < Ntime; i_t++) {

       // Lock all the coordinates except the passive linkage ones.
       for (const auto& coord : model.getComponentList<Coordinate>()) {

           auto state_traj = mt.getState(coord.getStateVariableNames()[0]);

           double v;
           if (IO::StartsWith(coord.getName(), "passive") ||
                   IO::StartsWith(coord.getName(), "motor")) {
               if (i_t == 0) {
                   v = coord.getDefaultValue();
               } else {
                   v = state_traj(i_t - 1);
               }
               coord.setLocked(s, false);
               coord.setValue(s, v, false);
               coord.setLocked(s, false);
           } else {
               v = state_traj(i_t);
               coord.setLocked(s, false);
               coord.setValue(s, v, false);
               coord.setLocked(s, true);
           }
       }

       // Resolve the kinematics.
       model.assemble(s);
       model.realizeVelocity(s);


       // Get and store the new linkage coordinates.
       for (const auto& coord : model.getComponentList<Coordinate>()) {
           if (IO::StartsWith(coord.getName(), "passive") ||
                   IO::StartsWith(coord.getName(), "motor")) {

               std::string val_name = coord.getStateVariableNames()[0];
               auto val = model.getStateVariableValue(s, val_name);
               auto val_traj = mt.getState(val_name);
               val_traj[i_t] = val;
               mt.setState(val_name, val_traj);

               std::string spd_name = coord.getStateVariableNames()[1];
               auto spd = model.getStateVariableValue(s, spd_name);
               auto spd_traj = mt.getState(spd_name);
               spd_traj[i_t] = spd;
               mt.setState(spd_name, spd_traj);
           }
       }
   }

   auto& time_vec = mt.getTime();
   // Initialize linkage coordinate velocities using forward differencing. 
   // Get and store the new linkage coordinates.
   for (const auto& coord : model.getComponentList<Coordinate>()) {
       if (IO::StartsWith(coord.getName(), "passive") ||
               IO::StartsWith(coord.getName(), "motor")) {

           std::string val_name = coord.getStateVariableNames()[0];
           auto val_traj = mt.getState(val_name);

           std::string spd_name = coord.getStateVariableNames()[1];
           auto spd_traj = mt.getState(spd_name);

           
           for (int i_t = 0; i_t < Ntime - 1; i_t++) {
               double dt = time_vec[i_t + 1] - time_vec[i_t];
               double spd = (val_traj[i_t + 1] - val_traj[i_t]) / dt;
               spd_traj[i_t] = spd;
           }

           // Backward Differencing.
           spd_traj[Ntime - 1] = (val_traj[Ntime - 1] - val_traj[Ntime - 2]) /
                                 (time_vec[Ntime - 1] - time_vec[Ntime - 2]);
           mt.setState(spd_name, spd_traj);
       }
   }

   return mt; 
}

// Set a gait prediction problem where the goal is to minimize effort (squared
// controls) over distance traveled while enforcing symmetry of the walking
// cycle and a prescribed average gait speed through endpoint constraints. The
// solution of the coordinate tracking problem is passed as an input argument
// and used as an initial guess for the prediction.
MocoSolution gaitPrediction(std::string model_file, MocoTrajectory guess,
        std::string track_file, double speed, double motor_weight,
        double track_weight, double speed_bound, double motor_bound,
        std::string scaling_method, int Nmesh, int Nparallel,
        std::string fn_prefix, int eff_exp_muscle, int eff_exp_motor, int NmaxIts,
        std::string motor_mode, double smooth, double k_pea,
        double t_pea, double gamma_lim, double lambda_lim) {

    using SimTK::Pi;

    MocoStudy study;
    study.setName(fn_prefix + "_gait_study");

    // Define the optimal control problem.
    // ===================================
    Model model_in(model_file);
    double l0, l1, l2, l3;
    l0 = 0.07;
    l1 = 0.05;
    l2 = 0.11;
    l3 = 0.1;
    configureLinkage(
            l0, l1, l2, l3, model_in, -Pi / 180.0 * 5.0, -Pi / 180.0 * 175.0);
    model_in.print("test.osim");

    bool has_spring = false;
    for (auto& spring : model_in.updComponentList<SpringGeneralizedForce>()) {
        std::cout << "Found Spring: " << spring.getName() << std::endl;
        has_spring = true;
        spring.setStiffness(k_pea);
        spring.setRestLength(t_pea);
    }

    MocoProblem& problem = study.updProblem();
    ModelProcessor modelprocessor = ModelProcessor(model_in);
    problem.setModelProcessor(modelprocessor);
    Model model = modelprocessor.process();
    model.initSystem();

    // Detect motors: Torque Actuator, DCMotor, or neither.
    bool has_dcmotor = false;
    for (const auto& motor : model.getComponentList<DCMotor>()) {
        std::cout << "Found DCMotor: " << motor.getName() << std::endl;
        has_dcmotor = true;
    }

    bool has_idealtorque = false;
    for (const auto& motor : model.getComponentList<TorqueActuator>()) {
        std::cout << "Found Torque Actuator: " << motor.getName() << std::endl;
        has_idealtorque = true;
    }

    

    // Goals.
    // =====
    // Symmetry.
    auto* symmetryGoal =
            problem.addGoal<MocoPeriodicityGoal>("symmetry_constraint");

    // Symmetric coordinate values (except for pelvis_tx) and speeds.
    for (const auto& coord : model.getComponentList<Coordinate>()) {
        if (IO::StartsWith(coord.getName(), "passive") ||
                IO::StartsWith(coord.getName(), "motor")) {
            // Redundant to make these symmteric.
            continue;
        } else if (IO::EndsWith(coord.getName(), "_r")) {
            symmetryGoal->addStatePair({coord.getStateVariableNames()[0],
                    std::regex_replace(coord.getStateVariableNames()[0],
                            std::regex("_r"), "_l")});
            symmetryGoal->addStatePair({coord.getStateVariableNames()[1],
                    std::regex_replace(coord.getStateVariableNames()[1],
                            std::regex("_r"), "_l")});
        } else if (IO::EndsWith(coord.getName(), "_l")) {
            symmetryGoal->addStatePair({coord.getStateVariableNames()[0],
                    std::regex_replace(coord.getStateVariableNames()[0],
                            std::regex("_l"), "_r")});
            symmetryGoal->addStatePair({coord.getStateVariableNames()[1],
                    std::regex_replace(coord.getStateVariableNames()[1],
                            std::regex("_l"), "_r")});
        } else if (!IO::EndsWith(coord.getName(), "_l") &&
                   !IO::EndsWith(coord.getName(), "_r") &&
                   !IO::EndsWith(coord.getName(), "_tx")) {
            symmetryGoal->addStatePair({coord.getStateVariableNames()[0],
                    coord.getStateVariableNames()[0]});
            symmetryGoal->addStatePair({coord.getStateVariableNames()[1],
                    coord.getStateVariableNames()[1]});
        }
    }
    symmetryGoal->addStatePair({"/jointset/groundPelvis/pelvis_tx/speed"});

    // Symmetric motor current if DC motor.
    for (const auto& motor : model.getComponentList<DCMotor>()) {
        std::string s_name = motor.getAbsolutePathString() + "/current";
        if (IO::EndsWith(motor.getName(), "_r")) {
            symmetryGoal->addStatePair({s_name,
                    std::regex_replace(s_name, std::regex("_r"), "_l")});
        } else if (IO::EndsWith(motor.getName(), "_l")) {
            symmetryGoal->addStatePair({s_name,
                    std::regex_replace(s_name, std::regex("_l"), "_r")});
        }
    }

    // Symmetric coordinate actuator controls.
    symmetryGoal->addControlPair({"/lumbarAct"});

    if (has_dcmotor || has_idealtorque) {
        symmetryGoal->addControlPair({"/motor_r", "/motor_l"});
        symmetryGoal->addControlPair({"/motor_l", "/motor_r"});
    }

    // Symmetric muscle activations.
    for (const auto& muscle : model.getComponentList<Muscle>()) {
        if (IO::EndsWith(muscle.getName(), "_r")) {
            symmetryGoal->addStatePair({muscle.getStateVariableNames()[0],
                    std::regex_replace(muscle.getStateVariableNames()[0],
                            std::regex("_r"), "_l")});
        }
        if (IO::EndsWith(muscle.getName(), "_l")) {
            symmetryGoal->addStatePair({muscle.getStateVariableNames()[0],
                    std::regex_replace(muscle.getStateVariableNames()[0],
                            std::regex("_l"), "_r")});
        }
    }
    //// Prescribed average gait speed.
    auto* speedGoal = problem.addGoal<MocoAverageSpeedGoal>("speed");
    speedGoal->set_desired_average_speed(speed);

    // Motor effort over time.
    if (has_dcmotor || has_idealtorque) {
        if (motor_mode == "exponent") {

            auto* motorEffortGoal =
                    problem.addGoal<MocoControlGoal>("effort_motor", 10);
            motorEffortGoal->setExponent(eff_exp_motor);
            motorEffortGoal->setDivideByDuration(true);
            motorEffortGoal->setWeightForControlPattern(
                    ".*", 0); // Set everything to 0 first.
            motorEffortGoal->setWeightForControl("/motor_r", motor_weight);
            motorEffortGoal->setWeightForControl("/motor_l", motor_weight);
        } else if (motor_mode == "energy") {

            auto* motorEnergyGoal = problem.addGoal<MocoEnergyGoal>(
                    "energy_motor", motor_weight);
            motorEnergyGoal->setDivideByDuration(true);

            if (has_idealtorque) {
                // Torque and speed.
                motorEnergyGoal->addPair(
                        {"/motor_r", "/j_in_r/motor_angle_r/speed", true});
                motorEnergyGoal->addPair(
                        {"/motor_l", "/j_in_l/motor_angle_l/speed", true});
            } else if (has_dcmotor) {
                // Voltage and current.
                motorEnergyGoal->addPair({"/motor_r", "/motor_r/current"});
                motorEnergyGoal->addPair({"/motor_l", "/motor_l/current"});
            }
            motorEnergyGoal->setSmoothScale(smooth);
        } else {
            OPENSIM_THROW(InvalidArgument,
                    "Invalid motor mode: '" + motor_mode +
                            "'. Must be 'exponent' or 'energy.'");
        }
    }

    // Effort over time.
    auto* muscleEffortGoal =
            problem.addGoal<MocoControlGoal>("effort_muscle", 10);
    muscleEffortGoal->setExponent(eff_exp_muscle);
    muscleEffortGoal->setDivideByDuration(true);

    if (has_dcmotor || has_idealtorque) {
        muscleEffortGoal->setWeightForControl("/motor_r", 0);
        muscleEffortGoal->setWeightForControl("/motor_l", 0);
    }

    // State Tracking (optional)
    if (track_file.length() > 0) {
        auto* stateTracking = problem.addGoal<MocoStateTrackingGoal>(
                "tracking", track_weight);
        stateTracking->setAllowUnusedReferences(true);
        auto ref_table_proc =
                TableProcessor(track_file) | TabOpLowPassFilter(20);
        auto ref_table = ref_table_proc.process(&model);
        stateTracking->setReference(ref_table);
        STOFileAdapter::write(ref_table, "track_ref.sto");
        auto track_traj = MocoTrajectory(track_file);
        problem.setTimeInfo(
                track_traj.getInitialTime(), track_traj.getFinalTime());
    } else {
        problem.setTimeInfo(0.0, {0.2, 0.7});
        // problem.setTimeInfo(0.0, 0.457201469);
        // problem.setTimeInfo(0.0, 0.4);
    }

    // Bounds.
    // =======
    problem.setStateInfo("/jointset/groundPelvis/pelvis_tilt/value",
            {-20 * Pi / 180, -10 * Pi / 180});
    problem.setStateInfo(
            "/jointset/groundPelvis/pelvis_tx/value", {0, 1}, 0, {0.2, 1.0});
    problem.setStateInfo(
            "/jointset/groundPelvis/pelvis_ty/value", {0.80, 1.00});
    problem.setStateInfo("/jointset/hip_l/hip_flexion_l/value",
            {-10 * Pi / 180, 60 * Pi / 180});
    problem.setStateInfo("/jointset/hip_r/hip_flexion_r/value",
            {-10 * Pi / 180, 60 * Pi / 180});
    problem.setStateInfo(
            "/jointset/knee_l/knee_angle_l/value", {-90 * Pi / 180, 0});
    problem.setStateInfo(
            "/jointset/knee_r/knee_angle_r/value", {-90 * Pi / 180, 0});
    //problem.setStateInfo("/jointset/ankle_l/ankle_angle_l/value",
    //        {-25 * Pi / 180, 25 * Pi / 180});
    //problem.setStateInfo("/jointset/ankle_r/ankle_angle_r/value",
    //        {-25 * Pi / 180, 25 * Pi / 180});
    problem.setStateInfo("/jointset/lumbar/lumbar/value", {0, 20 * Pi / 180});

    double s_mech = 200;
    double s_ankle = speed_bound;

    problem.setStateInfo(
            "/jointset/ankle_r/ankle_angle_r/speed", {-s_ankle, s_ankle});
    problem.setStateInfo(
            "/jointset/ankle_l/ankle_angle_l/speed", {-s_ankle, s_ankle});

    if (auto body = model.findComponent<Coordinate>("motor_angle")) {
        problem.setStateInfo("/j_in_r/motor_angle_r/speed", {-s_mech, s_mech});
        problem.setStateInfo("/j_in_l/motor_angle_l/speed", {-s_mech, s_mech});
    }

    if (auto body = model.findComponent<Coordinate>("passive_angle1")) {
        problem.setStateInfo(
                "/j_passive1_r/passive_angle1_r/speed", {-s_mech, s_mech});
        problem.setStateInfo(
                "/j_passive1_l/passive_angle1_l/speed", {-s_mech, s_mech});
    }
    problem.setStateInfo("/jointset/groundPelvis/pelvis_tx/speed", {-5, 5});
    problem.setStateInfo("/jointset/groundPelvis/pelvis_ty/speed", {-5, 5});

    if (auto body = model.findComponent<Coordinate>("SEA_angle")) {
        problem.setStateInfo("/j_in_r/SEA_angle_r/value", {-3, 3});
        problem.setStateInfo("/j_in_l/SEA_angle_l/value", {-3, 3});
    }

    // If the model uses DCMotor actuators, set their current limits.
    for (const auto& motor : model.getComponentList<DCMotor>()) {
        double max_current = motor.getMaximumCurrent();
        problem.setStateInfo(motor.getAbsolutePathString() + "/current",
                {-max_current, max_current});

        // Do not enforce current limits in dynamic equations.
        motor.get_enforce_current_limit(false);
    }

    problem.setMultiplierBounds({-lambda_lim, lambda_lim});

    if (has_dcmotor || has_idealtorque) {
        problem.setControlInfo(
                "/motor_r", {-motor_bound, motor_bound}, {}, {}, 1);
        problem.setControlInfo(
                "/motor_l", {-motor_bound, motor_bound}, {}, {}, 1);
    }
    problem.setControlInfo("/lumbarAct", {-0.5, 0.5}, {}, {}, 1);

    // problem.setMultiplierScaler(1);

    // Configure the solver.
    // =====================
    auto& solver = study.initCasADiSolver();
    solver.set_num_mesh_intervals(Nmesh);
    solver.set_verbosity(2);
    solver.set_optim_solver("ipopt");
    solver.set_optim_convergence_tolerance(1e-4);
    solver.set_optim_constraint_tolerance(1e-4);
    solver.set_optim_max_iterations(NmaxIts);
    solver.set_enforce_constraint_derivatives(true);
    solver.set_minimize_lagrange_multipliers(false);
    solver.set_velocity_correction_bounds({-gamma_lim, gamma_lim});
    solver.set_optim_ipopt_opt_filename(fn_prefix + ".opt");

    // Set the guess.
    // ==============
    auto processed_guess = genInitialGuess(study, guess);
    if (track_file.length() > 0) {
        auto &linkage_processed_guess =
                generateLinkageAngles(model, processed_guess);
        solver.setGuess(linkage_processed_guess);
    } else {
        solver.setGuess(processed_guess);
    }
    solver.set_scaling_method(scaling_method);
    solver.set_parallel(Nparallel);
    solver.set_optim_finite_difference_scheme("central");
    solver.getGuess().write(fn_prefix + "_guess.sto");

    // Solve problem.
    // ==============
    study.print(fn_prefix + "_study.moco");
    MocoSolution solution = study.solve();
    solution.unseal();
    solution.write(fn_prefix + "_solution_halfcycle.sto");
    auto full = createPeriodicTrajectory(solution);
    full.write(fn_prefix + "_solution_fullcycle.sto");

    // Extract ground reaction forces.
    // ===============================
    std::vector<std::string> contact_all;
    std::vector<std::string> contact_r;
    std::vector<std::string> contact_l;
    contact_r.push_back("contactHeel_r");
    contact_r.push_back("contactFront_r");
    contact_l.push_back("contactHeel_l");
    contact_l.push_back("contactFront_l");
    contact_all = contact_r;
    contact_all.insert(std::end(contact_all), std::begin(contact_l), std::end(contact_l));
    TimeSeriesTable externalForcesTableFlat =
            createExternalLoadsTableForGait(model, full, contact_r, contact_l);
    STOFileAdapter::write(
            externalForcesTableFlat, fn_prefix +"_solutionGRF_fullcycle.sto");

    TimeSeriesTable contactForcesTableFlat =
            createContactForceTableForGait(model, full, contact_all);
    STOFileAdapter::write(
            contactForcesTableFlat, fn_prefix +"_solutionContact_fullcycle.sto");

    // study.visualize(full);

    return solution;
}

int main(int argc, char* argv[]) {
    try {
        std::string model_file, guess_file, track_file, scaling_method,
                motor_mode;
        double motor_weight, track_weight, speed_bound, motor_bound, speed,
                smooth, k_pea, t_pea, gamma_lim, lambda_lim;
        int Nmesh, Nparallel, eff_exp_motor, eff_exp_muscle, NmaxIts;

        // Log the version for reference.
        std::cout << "2D Gait Generation -- Tyler Morrison 2021" << std::endl;
        std::cout << "Version compiled " << __DATE__ << " at " __TIME__ << "."
                  << std::endl;
        bool predict_only;
        if (argc == 21) {
            std::cout << "Run in two-step track + predict mode..." << std::endl;
            predict_only = false;
            model_file = argv[1];
            guess_file = argv[2];
            track_file = argv[3];
            motor_weight = atof(argv[4]);
            track_weight = atof(argv[5]);
            speed_bound = atof(argv[6]);
            motor_bound = atof(argv[7]);
            scaling_method = argv[8];
            Nmesh = atoi(argv[9]);
            Nparallel = atoi(argv[10]);
            speed = atof(argv[11]);
            eff_exp_muscle = atoi(argv[12]);
            eff_exp_motor = atoi(argv[13]);
            NmaxIts = atoi(argv[14]);
            motor_mode = argv[15];
            smooth = atof(argv[16]);
            k_pea = atof(argv[17]);
            t_pea = atof(argv[18]);
            gamma_lim = atof(argv[19]);
            lambda_lim = atof(argv[20]);
        } else if (argc == 19) {
            std::cout << "Run in one-step predict mode..." << std::endl;
            predict_only = true;
            model_file = argv[1];
            guess_file = argv[2];
            track_file = "";
            motor_weight = atof(argv[3]);
            track_weight = 0.0;
            speed_bound = atof(argv[4]);
            motor_bound = atof(argv[5]);
            scaling_method = argv[6];
            Nmesh = atoi(argv[7]);
            Nparallel = atoi(argv[8]);
            speed = atof(argv[9]);
            eff_exp_muscle = atoi(argv[10]);
            eff_exp_motor = atoi(argv[11]);
            NmaxIts = atoi(argv[12]);
            motor_mode = argv[13];
            smooth = atof(argv[14]);
            k_pea = atof(argv[15]);
            t_pea = atof(argv[16]);
            gamma_lim = atof(argv[17]);
            lambda_lim = atof(argv[18]);
        } else if (argc == 3) {
            std::cout << "Generate contact force mode..." << std::endl;
            model_file = argv[1];
            std::string trajectory_file = argv[2];

            MocoTrajectory traj(trajectory_file);
            MocoStudy study;
            study.setName("visualizer");
            Model model_in(model_file);
            MocoProblem& problem = study.updProblem();
            problem.setModelProcessor(ModelProcessor(model_in));


            // Extract ground reaction forces.
            // ===============================
            std::vector<std::string> contact_all;
            contact_all.push_back("contactHeel_r");
            contact_all.push_back("contactFront_r");
            contact_all.push_back("contactHeel_l");
            contact_all.push_back("contactFront_l");

            auto split_pos = trajectory_file.find("_");
            std::string prefix = trajectory_file.substr(0, split_pos);

            TimeSeriesTable contactForcesTableFlat =
                    createContactForceTableForGait(
                            model_in, traj, contact_all);

            std::string fn = prefix + "_solutionContact_fullcycle.sto";
            std::cout << "Writing contact forces to " << fn << std::endl;
            STOFileAdapter::write(contactForcesTableFlat, fn);

            return EXIT_SUCCESS;

        } else {
            std::cerr << "Input parse failure." << std::endl;
            return EXIT_FAILURE;
        }

        std::cout << "***************************************" << std::endl;
        std::cout << "1: OSIM Model: " << model_file << std::endl;
        std::cout << "2: Guess STO file: " << guess_file << std::endl;
        std::cout << "3: Track STO file: " << track_file << std::endl;
        std::cout << "4: Motor Weight: " << motor_weight << std::endl;
        std::cout << "5: Track Weight: " << track_weight << std::endl;
        std::cout << "6: Ankle speed bound: " << speed_bound << std::endl;
        std::cout << "7: Motor control bound: " << motor_bound << std::endl;
        std::cout << "8: Scaling methods: " << scaling_method << std::endl;
        std::cout << "9: Nmesh: " << Nmesh << std::endl;
        std::cout << "10: Nparallel: " << Nparallel << std::endl;
        std::cout << "11: Walking Speed: " << speed << std::endl;
        std::cout << "12: Muscle Effort Exponent: " << eff_exp_muscle << std::endl;
        std::cout << "13: Motor Effort Exponent: " << eff_exp_motor << std::endl;
        std::cout << "14: Max Iterations: " << NmaxIts << std::endl;
        std::cout << "15: Motor Cost Mode: " << motor_mode << std::endl;
        std::cout << "16: Energy smoothing constant: " << smooth << std::endl;
        std::cout << "17: Parallel spring constant: " << k_pea << std::endl;
        std::cout << "18: Parallel spring equilibrium: " << t_pea << std::endl;
        std::cout << "19: Velocity correction limit: " << gamma_lim << std::endl;
        std::cout << "20: Constraint force limit: " << lambda_lim << std::endl;
        std::cout << "***************************************" << std::endl;

        // Solve the initial tracking problem.
        MocoTrajectory guess_traj(guess_file);

        if (predict_only) {
            gaitPrediction(model_file, guess_traj, "", speed, motor_weight, 0,
                    speed_bound, motor_bound, scaling_method, Nmesh, Nparallel,
                    "predict", eff_exp_muscle, eff_exp_motor, NmaxIts,
                    motor_mode, smooth, k_pea, t_pea, gamma_lim, lambda_lim);
        } else {
        auto track_sol = gaitPrediction(model_file, guess_traj, track_file,
                speed, motor_weight, track_weight, speed_bound, motor_bound,
                scaling_method, Nmesh, Nparallel, "track", eff_exp_muscle, eff_exp_motor,
                NmaxIts, motor_mode, smooth, k_pea, t_pea, gamma_lim, lambda_lim);

        // Solve the pure prediction problem.
        // Ignores the track file and the track weight.
        gaitPrediction(model_file, track_sol, "", speed, motor_weight, 0,
                speed_bound, motor_bound, scaling_method, Nmesh, Nparallel,
                "predict", eff_exp_muscle, eff_exp_motor, NmaxIts, motor_mode,
                smooth, k_pea, t_pea, gamma_lim, lambda_lim);
        }

    } catch (const std::exception& e) { std::cout << e.what() << std::endl; }
    return EXIT_SUCCESS;
}
