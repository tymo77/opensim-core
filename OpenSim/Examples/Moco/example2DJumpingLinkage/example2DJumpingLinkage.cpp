/* -------------------------------------------------------------------------- *
 * OpenSim Moco: example2DJumping.cpp                                         *
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

//std::string SCALING_METHOD = "none";
//std::string SCALING_METHOD = "info";
std::string SCALING_METHOD = "bounds";

using namespace OpenSim;

void jumpPrediction(double motor_weight, double motor_scaler, double slack_scaler, double mult_scaler) {

    using SimTK::Pi;

    MocoStudy study;
    study.setName("jumpPrediction");

    // Define the optimal control problem.
    // ===================================
    MocoProblem& problem = study.updProblem();
    ModelProcessor modelprocessor = ModelProcessor("2D_linkage_leg.osim");
    problem.setModelProcessor(modelprocessor);
     Model model = modelprocessor.process();
     model.initSystem();

    // Goals.
    // =====

    // Effort over time.
    auto* motorEffortGoal = problem.addGoal<MocoControlGoal>("effort_motor", 10);
    motorEffortGoal->setExponent(2);
    motorEffortGoal->setDivideByDuration(true);
    motorEffortGoal->setWeightForControlPattern(".*", 0); // Set everything to 0 first.
    motorEffortGoal->setWeightForControl("/motor_r", motor_weight);

    // Effort over time.
    auto* muscleEffortGoal = problem.addGoal<MocoControlGoal>("effort_muscle", 10);
    muscleEffortGoal->setExponent(2);
    muscleEffortGoal->setDivideByDuration(true);
    muscleEffortGoal->setWeightForControl("/motor_r", 0);

    // Bounds.
    // =======
    problem.setTimeInfo(0, {0.1, 1.5});
    problem.setStateInfo(
            "/jointset/groundPelvis/pelvis_ty/value", {0.50, 1.30}, 0.88, 1.25);
    problem.setStateInfo("/jointset/hip_r/hip_flexion_r/value",
            {-10.0 * Pi / 180.0, 60.0 * Pi / 180.0}, 25.0 * Pi / 180.0);
    problem.setStateInfo("/jointset/knee_r/knee_angle_r/value",
            {-90.0 * Pi / 180, 0.0}, -50.0 * Pi / 180.0);
    problem.setStateInfo("/jointset/ankle_r/ankle_angle_r/value",
            {-25.0 * Pi / 180.0, 25.0 * Pi / 180.0}, 25.0 * Pi / 180.0);

    problem.setStateInfo("/jointset/ankle_r/ankle_angle_r/speed", {-200, 200});

    problem.setStateInfo("/j_in_r/motor_angle_r/speed", {-200, 200});

    problem.setStateInfo("/j_passive1_r/passive_angle1_r/speed", {-200, 200});

    problem.setStateInfo("/jointset/groundPelvis/pelvis_ty/speed", {-10, 10}, 0, 0);


    problem.setMultiplierBounds({-5000, 5000});
    problem.setControlInfo("/motor_r", {-200, 200}, {}, {}, motor_scaler);
    
    problem.setMultiplierScaler(mult_scaler);
    

    // Configure the solver.
    // =====================
    auto& solver = study.initCasADiSolver();
    solver.set_num_mesh_intervals(30);
    solver.set_verbosity(2);
    solver.set_optim_solver("ipopt");
    solver.set_optim_convergence_tolerance(1e-4);
    solver.set_optim_constraint_tolerance(1e-4);
    solver.set_optim_max_iterations(10000);
    solver.set_enforce_constraint_derivatives(true);
    solver.set_minimize_lagrange_multipliers(false);
    //solver.set_lagrange_multiplier_weight(0.1);
    solver.set_velocity_correction_bounds({-100, 100});
    // Use the solution from the tracking simulation as initial guess.
    solver.set_scaling_method(SCALING_METHOD);
    solver.set_velocity_correction_scaler(slack_scaler);
    solver.set_parallel(1);
    //solver.setGuess(solver.createGuess("bounds"));
    solver.setGuessFile("prev_sol.sto");
    solver.getGuess().write("guess.sto");
    
    // Solve problem.
    // ==============
    MocoSolution solution = study.solve();
    solution.unseal();
    solution.write("jumpPrediction.sto");

    // Extract ground reaction forces.
    // ===============================
    std::vector<std::string> contact_r, contact_l;
    contact_r.push_back("contactHeel_r");
    contact_r.push_back("contactFront_r");
    TimeSeriesTable externalForcesTableFlat = createExternalLoadsTableForGait(
            model, solution, contact_r, contact_l);
    STOFileAdapter::write(externalForcesTableFlat,
            "jumpPrediction_GRF.sto");

    study.visualize(solution);
}

int main(int argc, char* argv[]) {
    try {

        /*double motor_weight = atof(argv[1]);
        double motor_scaler = atof(argv[2]);
        double slack_scaler = atof(argv[3]);
        double mult_scaler = atof(argv[4]);

        std::cout << "Motor Weight: " << motor_weight << std::endl;
        std::cout << "Motor Scaler: " << motor_scaler << std::endl;
        std::cout << "Slack Scaler: " << slack_scaler << std::endl;
        std::cout << "Multiplier Scaler: " << mult_scaler << std::endl;*/

        //const MocoSolution gaitTrackingSolution = gaitTracking();
        //gaitPrediction(motor_weight, motor_scaler, slack_scaler, mult_scaler);
        jumpPrediction(1e-3, 1, 1, 1);
    } catch (const std::exception& e) { std::cout << e.what() << std::endl; }
    return EXIT_SUCCESS;
}
