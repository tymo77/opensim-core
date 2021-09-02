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

#include <OpenSim/Common/STOFileAdapter.h>
#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/Common/TableUtilities.h>
#include <OpenSim/Actuators/DCMotor.h>
#include <OpenSim/Actuators/TorqueActuator.h>

using namespace OpenSim;

// Set a gait prediction problem where the goal is to minimize effort (squared
// controls) over distance traveled while enforcing symmetry of the walking
// cycle and a prescribed average gait speed through endpoint constraints. The
// solution of the coordinate tracking problem is passed as an input argument
// and used as an initial guess for the prediction.
MocoSolution gaitPrediction(std::string model_file, MocoTrajectory guess,
        std::string track_file, double speed, double motor_weight, double track_weight,
        double speed_bound, double motor_bound,
        std::string scaling_method, int Nmesh, int Nparallel,
        std::string fn_prefix, int eff_exp, int NmaxIts, std::string motor_mode, double smooth) {

    using SimTK::Pi;

    MocoStudy study;
    study.setName(fn_prefix + "_gait_study");

    // Define the optimal control problem.
    // ===================================
    Model model_in(model_file);
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
    auto* symmetryGoal = problem.addGoal<MocoPeriodicityGoal>("symmetry_constraint");
    
    // Symmetric coordinate values (except for pelvis_tx) and speeds.
    for (const auto& coord : model.getComponentList<Coordinate>()) {
        if (IO::StartsWith(coord.getName(), "passive") ||
                IO::StartsWith(coord.getName(), "motor")) {
            // Redundant to make these symmteric.
            continue;
        }
        else if (IO::EndsWith(coord.getName(), "_r")) {
            symmetryGoal->addStatePair({coord.getStateVariableNames()[0],
                    std::regex_replace(coord.getStateVariableNames()[0],
                            std::regex("_r"), "_l")});
            symmetryGoal->addStatePair({coord.getStateVariableNames()[1],
                    std::regex_replace(coord.getStateVariableNames()[1],
                            std::regex("_r"), "_l")});
        }
        else if (IO::EndsWith(coord.getName(), "_l")) {
            symmetryGoal->addStatePair({coord.getStateVariableNames()[0],
                    std::regex_replace(coord.getStateVariableNames()[0],
                            std::regex("_l"), "_r")});
            symmetryGoal->addStatePair({coord.getStateVariableNames()[1],
                    std::regex_replace(coord.getStateVariableNames()[1],
                            std::regex("_l"), "_r")});
        }
        else if (!IO::EndsWith(coord.getName(), "_l") &&
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
    // Prescribed average gait speed.
    auto* speedGoal = problem.addGoal<MocoAverageSpeedGoal>("speed");
    speedGoal->set_desired_average_speed(speed);

    // Motor effort over time.
    if (has_dcmotor || has_idealtorque) {
        if (motor_mode == "exponent") {

            auto* motorEffortGoal =
                    problem.addGoal<MocoControlGoal>("effort_motor", 10);
            motorEffortGoal->setExponent(eff_exp);
            motorEffortGoal->setDivideByDuration(true);
            motorEffortGoal->setWeightForControlPattern(
                    ".*", 0); // Set everything to 0 first.
            motorEffortGoal->setWeightForControl("/motor_r", motor_weight);
            motorEffortGoal->setWeightForControl("/motor_l", motor_weight);
        } else if (motor_mode == "energy") {

            auto* motorEnergyGoal = problem.addGoal<MocoEnergyGoal>(
                    "energy_motor", motor_weight);
            motorEnergyGoal->setDivideByDuration(true);
            motorEnergyGoal->addPair({"/motor_r", "/motor_r/current"});
            motorEnergyGoal->addPair({"/motor_l", "/motor_l/current"});
            motorEnergyGoal->setSmoothScale(smooth);
        } else {
            OPENSIM_THROW(InvalidArgument,
                    "Invalid motor mode: '" + motor_mode +
                            "'. Must be 'exponent' or 'energy.'");
        }
    }

    // Effort over time.
    auto* muscleEffortGoal = problem.addGoal<MocoControlGoal>("effort_muscle", 10);
    muscleEffortGoal->setExponent(eff_exp);
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
        auto ref_table_proc = TableProcessor(track_file) | TabOpLowPassFilter(20);
        auto ref_table = ref_table_proc.process(&model);
        stateTracking->setReference(ref_table);
        STOFileAdapter::write(ref_table, "track_ref.sto");
        auto track_traj = MocoTrajectory(track_file);
        problem.setTimeInfo(
                track_traj.getInitialTime(), track_traj.getFinalTime());
    } else {
        problem.setTimeInfo(0.0, {0.2, 0.7});
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
    problem.setStateInfo("/jointset/ankle_l/ankle_angle_l/value",
            {-25 * Pi / 180, 25 * Pi / 180});
    problem.setStateInfo("/jointset/ankle_r/ankle_angle_r/value",
            {-25 * Pi / 180, 25 * Pi / 180});
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


    // If the model uses DCMotor actuators, set their current limits.
    for (const auto& motor : model.getComponentList<DCMotor>()) {
        double max_current = motor.getMaximumCurrent();
        problem.setStateInfo(motor.getAbsolutePathString()  + "/current",
                {-max_current, max_current});
        
        //Do not enforce current limits in dynamic equations.
        motor.get_enforce_current_limit(false);
    }

    problem.setMultiplierBounds({-3000, 3000});

    if (has_dcmotor || has_idealtorque) {
        problem.setControlInfo(
                "/motor_r", {-motor_bound, motor_bound}, {}, {}, 1);
        problem.setControlInfo(
                "/motor_l", {-motor_bound, motor_bound}, {}, {}, 1);
    }
    problem.setControlInfo("/lumbarAct", {-0.5, 0.5}, {}, {}, 1);
    
    //problem.setMultiplierScaler(1);
    

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
    solver.set_velocity_correction_bounds({-1, 1});
    solver.set_optim_ipopt_opt_filename(fn_prefix + ".opt");
    // Use the solution from the tracking simulation as initial guess.
    solver.setGuess(guess);
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
    std::vector<std::string> contact_r;
    std::vector<std::string> contact_l;
    contact_r.push_back("contactHeel_r");
    contact_r.push_back("contactFront_r");
    contact_l.push_back("contactHeel_l");
    contact_l.push_back("contactFront_l");
    TimeSeriesTable externalForcesTableFlat = createExternalLoadsTableForGait(
            model, full, contact_r, contact_l);
    STOFileAdapter::write(externalForcesTableFlat,
        fn_prefix +  + "_solutionGRF_fullcycle.sto");

    //study.visualize(full);

    return solution;
}

int main(int argc, char* argv[]) {
    try {
        std::string model_file, guess_file, track_file, scaling_method, motor_mode;
        double motor_weight, track_weight, speed_bound, motor_bound, speed, smooth;
        int Nmesh, Nparallel, eff_exp, NmaxIts;

        // Log the version for reference.
        std::cout << "2D Gait Generation -- Tyler Morrison 2021"
                  << std::endl;
        std::cout << "Version compiled " << __DATE__ << " at " __TIME__ << "."
                  << std::endl;

        if (argc == 16) {
            std::cout << "Run in two-step track + predict mode..." << std::endl;
            model_file       = argv[1];
            guess_file       = argv[2];
            track_file       = argv[3];
            motor_weight     = atof(argv[4]);
            track_weight     = atof(argv[5]);
            speed_bound      = atof(argv[6]);
            motor_bound     = atof(argv[7]);
            scaling_method   = argv[8];
            Nmesh            = atoi(argv[9]);
            Nparallel        = atoi(argv[10]);
            speed            = atof(argv[11]);
            eff_exp          = atoi(argv[12]);
            NmaxIts          = atoi(argv[13]);
            motor_mode       = argv[14];
            smooth           = atof(argv[15]);
        } else if (argc == 14) {
            std::cout << "Run in one-step predict mode..." << std::endl;
            model_file       = argv[1];
            guess_file       = argv[2];
            track_file       = "";
            motor_weight     = atof(argv[3]);
            track_weight     = 0.0;
            speed_bound      = atof(argv[4]);
            motor_bound     = atof(argv[5]);
            scaling_method   = argv[6];
            Nmesh            = atoi(argv[7]);
            Nparallel        = atoi(argv[8]);
            speed            = atof(argv[9]);
            eff_exp          = atoi(argv[10]);
            NmaxIts          = atoi(argv[11]);
            motor_mode       = argv[12];
            smooth           = atof(argv[13]);
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
        std::cout << "12: Effort Exponent: " << eff_exp << std::endl;
        std::cout << "13: Max Iterations: " << NmaxIts << std::endl;
        std::cout << "14: Motor Cost Mode: " << motor_mode << std::endl;
        std::cout << "15: Energy smoothing constant: " << smooth << std::endl;
        std::cout << "***************************************" << std::endl;

        // Solve the initial tracking problem.
        MocoTrajectory guess_traj(guess_file);

        auto track_sol = gaitPrediction(
            model_file, guess_traj, track_file,
            speed, motor_weight, track_weight, speed_bound, motor_bound, scaling_method, Nmesh, Nparallel, "track",
                eff_exp, NmaxIts, motor_mode, smooth);

        // Solve the pure prediction problem.
        // Ignores the track file and the track weight.
        gaitPrediction(model_file, track_sol, "",
            speed, motor_weight, 0, speed_bound, motor_bound, scaling_method, Nmesh, Nparallel,
                "predict", eff_exp, NmaxIts, motor_mode, smooth);

    } catch (const std::exception& e) { std::cout << e.what() << std::endl; }
    return EXIT_SUCCESS;
}
