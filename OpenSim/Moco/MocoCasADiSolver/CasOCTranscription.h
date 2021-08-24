#ifndef OPENSIM_CASOCTRANSCRIPTION_H
#define OPENSIM_CASOCTRANSCRIPTION_H
/* -------------------------------------------------------------------------- *
 * OpenSim: CasOCTranscription.h                                              *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2018 Stanford University and the Authors                     *
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

#include "CasOCSolver.h"

namespace CasOC {

/// This is the base class for transcription schemes that convert a
/// CasOC::Problem into a general nonlinear programming problem. If you are
/// creating a new derived class, make sure to override all virtual functions
/// and obey the settings that the user specified in the CasOC::Solver.
class Transcription {
public:
    Transcription(const Solver& solver, const Problem& problem)
            : m_solver(solver), m_problem(problem) {}
    virtual ~Transcription() = default;
    Iterate createInitialGuessFromBounds() const;
    /// Use the provided random number generator to generate an iterate.
    /// Random::Uniform is used if a generator is not provided. The generator
    /// should produce numbers with [-1, 1].
    Iterate createRandomIterateWithinBounds(
            const SimTK::Random* = nullptr) const;
    template <typename T>
    T createTimes(const T& initialTime, const T& finalTime) const {
        return (finalTime - initialTime) * m_grid + initialTime;
    }
    casadi::DM createQuadratureCoefficients() const {
        return createQuadratureCoefficientsImpl();
    }
    casadi::DM createMeshIndices() const {
        casadi::DM meshIndices = createMeshIndicesImpl();
        const auto shape = meshIndices.size();
        OPENSIM_THROW_IF(shape.first != 1 || shape.second != m_numGridPoints,
                OpenSim::Exception,
                "createMeshIndicesImpl() must return a row vector of shape "
                "length [1, {}], but a matrix of shape [{}, {}] was returned.",
                m_numGridPoints, shape.first, shape.second);
        OPENSIM_THROW_IF(!SimTK::isNumericallyEqual(
                                 casadi::DM::sum2(meshIndices).scalar(),
                                 m_numMeshPoints),
                OpenSim::Exception,
                "Internal error: sum of mesh indices should be the number of "
                "mesh points.");

        return meshIndices;
    }

    Solution solve(const Iterate& guessOrig);

    // These are the upper and lower bounds that variables send to the NLP
    // solver will be scaled between (if they have LB and UB).
    static double SCALED_LB;
    static double SCALED_UB;

protected:
    /// This must be called in the constructor of derived classes so that
    /// overridden virtual methods are accessible to the base class. This
    /// implementation allows initialization to occur during construction,
    /// avoiding an extra call on the instantiated object.
    /// pointsForInterpControls are grid points at which the transcription
    /// scheme applies constraints between control points.
    void createVariablesAndSetBounds(const casadi::DM& grid,
            int numDefectsPerMeshInterval,
            const casadi::DM& pointsForInterpControls = casadi::DM());

    /// We assume all functions depend on time and parameters.
    /// "inputs" is prepended by time and postpended (?) by parameters.
    casadi::MXVector evalOnTrajectory(const casadi::Function& pointFunction,
            const std::vector<Var>& inputs,
            const casadi::Matrix<casadi_int>& timeIndices) const;

    template <typename TRow, typename TColumn>
    void setVariableBounds(Var var, const TRow& rowIndices,
            const TColumn& columnIndices, const Bounds& bounds) {
        if (bounds.isSet()) {
            const auto& lower = bounds.lower;
            m_lowerBounds[var](rowIndices, columnIndices) = lower;
            const auto& upper = bounds.upper;
            m_upperBounds[var](rowIndices, columnIndices) = upper;
        } else {
            const auto inf = std::numeric_limits<double>::infinity();
            m_lowerBounds[var](rowIndices, columnIndices) = -inf;
            m_upperBounds[var](rowIndices, columnIndices) = inf;
        }
    }

    // Scale the given variables to the same range. Scale function
    // is a linear scale between the bounds. Bounds so close as to be equality
    // constraints are not scaled.
    template <typename TRow, typename TColumn>
    void setVariableScaleByBounds(Var var, const TRow& rowIndices,
            const TColumn& columnIndices, const Bounds& bounds) {
        if (bounds.isSet()) {
            const auto& lower = bounds.lower;
            const auto& upper = bounds.upper;

            // TODO: Best way to check for equality bounds?
            // Check if bounds are so close they are an eq. constraint.
            double m, b;
            if (abs(upper - lower) / abs(upper) > 1e-10 ||
                    abs(upper - lower) / abs(lower) > 1e-10 ||
                    abs(upper - lower) > 1e-10
                ) {
                m = (SCALED_UB - SCALED_LB) / (upper - lower);
                b = SCALED_LB - m * lower;
            } else {
                m = 1.0;
                b = (SCALED_LB + SCALED_UB)/2 - m * lower;
            }

            m_scalingMultiplier[var](rowIndices, columnIndices) = m;
            m_scalingShifter[var](rowIndices, columnIndices) = b;

        } else {
            m_scalingMultiplier[var](rowIndices, columnIndices) = 1.0;
            m_scalingShifter[var](rowIndices, columnIndices) = 0.0;
        }
    }

    template <typename TRow, typename TColumn>
    void setVariableScale(Var var, const TRow& rowIndices,
            const TColumn& columnIndices, const CasOCVarInfo& info) {

        ScaleMethod method = m_problem.getScaleMethod();
        switch (method) {
        case ScaleMethod::scale_off: 
            m_scalingMultiplier[var](rowIndices, columnIndices) = 1.0;
            m_scalingShifter[var](rowIndices, columnIndices) = 0.0;
            break;
        case ScaleMethod::scale_bounds: 
            setVariableScaleByBounds(var, rowIndices, columnIndices, info.bounds);
            break;
        case ScaleMethod::scale_var_info: 
            // TODO Is this a good test? Or is it too arbitrary?
            OPENSIM_THROW_IF(abs(info.variable_scaler) < 1e-10, OpenSim::Exception,
                "Scaling variable is too close to zero.")

            m_scalingMultiplier[var](rowIndices, columnIndices) = info.variable_scaler;
            m_scalingShifter[var](rowIndices, columnIndices) = 0.0;
            break;
        }

    }


    template <typename T>
    struct Constraints {
        T defects;
        T multibody_residuals;
        T auxiliary_residuals;
        T kinematic;
        std::vector<T> endpoint;
        std::vector<T> path;
        T interp_controls;
    };
    void printConstraintValues(const Iterate& it,
            const Constraints<casadi::DM>& constraints,
            std::ostream& stream = std::cout) const;
    void printObjectiveBreakdown(const Iterate& it,
            const casadi::DM& objectiveTerms,
            std::ostream& stream = std::cout) const;

    const Solver& m_solver;
    const Problem& m_problem;
    int m_numGridPoints = -1;
    int m_numMeshPoints = -1;
    int m_numMeshIntervals = -1;
    int m_numMeshInteriorPoints = -1;
    int m_numDefectsPerMeshInterval = -1;
    int m_numMultibodyResiduals = -1;
    int m_numAuxiliaryResiduals = -1;
    int m_numConstraints = -1;
    casadi::DM m_grid;
    casadi::DM m_pointsForInterpControls;
    casadi::MX m_times;
    casadi::MX m_duration;

private:
    VariablesMX m_vars;
    casadi::MX m_paramsTrajGrid;
    casadi::MX m_paramsTrajMesh;
    casadi::MX m_paramsTrajMeshInterior;
    VariablesDM m_lowerBounds;
    VariablesDM m_upperBounds;
    VariablesDM m_scalingMultiplier; // The multiplier used to scale the variable magnitude.
    VariablesDM m_scalingShifter; // The shift used to center the variable value
                                  // range.

    casadi::DM m_meshIndicesMap;
    casadi::Matrix<casadi_int> m_gridIndices;
    casadi::Matrix<casadi_int> m_meshIndices;
    casadi::Matrix<casadi_int> m_meshInteriorIndices;

    casadi::MX m_xdot; // State derivatives.

    casadi::MX m_objectiveTerms;
    std::vector<std::string> m_objectiveTermNames;

    Constraints<casadi::MX> m_constraints;
    Constraints<casadi::DM> m_constraintsLowerBounds;
    Constraints<casadi::DM> m_constraintsUpperBounds;

private:
    /// Override this function in your derived class to compute a vector of
    /// quadrature coeffecients (of length m_numGridPoints) required to set the
    /// the integral cost within transcribe().
    virtual casadi::DM createQuadratureCoefficientsImpl() const = 0;
    /// Override this function to specify the indicies in the grid where the
    /// mesh (or "knot") points lie.
    /// @note The returned vector must be a row vector of length m_numGridPoints
    /// with nonzero values at the mesh indices.
    virtual casadi::DM createMeshIndicesImpl() const = 0;
    /// Override this function in your derived class set the defect, kinematic,
    /// and path constraint errors required for your transcription scheme.
    virtual void calcDefectsImpl(const casadi::MX& x, const casadi::MX& xdot,
            casadi::MX& defects) const = 0;
    virtual void calcInterpolatingControlsImpl(const casadi::MX& /*controls*/,
            casadi::MX& /*interpControls*/) const {
        OPENSIM_THROW_IF(m_pointsForInterpControls.numel(), OpenSim::Exception,
                "Must provide constraints for interpolating controls.")
    }

    void transcribe();
    void setObjectiveAndEndpointConstraints();
    void calcDefects() {
        calcDefectsImpl(m_vars.at(states), m_xdot, m_constraints.defects);
    }
    void calcInterpolatingControls() {
        calcInterpolatingControlsImpl(
                m_vars.at(controls), m_constraints.interp_controls);
    }

    /// Use this function to ensure you iterate through variables in the same
    /// order.
    template <typename T>
    static std::vector<Var> getSortedVarKeys(const Variables<T>& vars) {
        std::vector<Var> keys;
        for (const auto& kv : vars) { keys.push_back(kv.first); }
        std::sort(keys.begin(), keys.end());
        return keys;
    }
    /// Convert the map of variables into a column vector, for passing onto
    /// nlpsol(), etc.
    template <typename T>
    static T flattenVariables(const CasOC::Variables<T>& vars) {
        std::vector<T> stdvec;
        for (const auto& key : getSortedVarKeys(vars)) {
            stdvec.push_back(vars.at(key));
        }
        return T::veccat(stdvec);
    }

    // Take the variables and linearly scale them to the fixed range.
    template <typename T>
    CasOC::Variables<T> scaleVariables(CasOC::Variables<T> vars) const {
        for (const auto& key : getSortedVarKeys(vars)) {
            vars.at(key) = scaleVar(key, vars.at(key));
        }
        return vars;
    }

    // Transform the variables back to their original values.
    template <typename T>
    CasOC::Variables<T> unScaleVariables(CasOC::Variables<T> vars) const {
        for (const auto& key : getSortedVarKeys(vars)) {
            vars.at(key) = unScaleVar(key, vars.at(key));
        }
        return vars;
    }

    CasOC::VariablesDM unScaleVariables(const CasOC::VariablesDM& vars) const {
        CasOC::VariablesDM out;
        for (const auto& key : getSortedVarKeys(vars)) {
            out[key] = unScaleVar(key, vars.at(key));
        }
        return out;
    }

    // Scale the variable by its linear bounds scaler.
    template <typename T>
    T scaleVar(CasOC::Var v, T x) const {
        return m_scalingMultiplier.at(v) * x + m_scalingShifter.at(v);
    }

    // Unscale a variable by its linear bounds scaler.
    template <typename T>
    T unScaleVar(CasOC::Var v, T x) const {
        return (x - m_scalingShifter.at(v)) / m_scalingMultiplier.at(v);
    }

    /// Convert the 'x' column vector into separate variables.
    CasOC::VariablesDM expandVariables(const casadi::DM& x) const {
        CasOC::VariablesDM out;
        using casadi::Slice;
        casadi_int offset = 0;
        for (const auto& key : getSortedVarKeys(m_vars)) {
            const auto& value = m_vars.at(key);
            // Convert a portion of the column vector into a matrix.
            out[key] = casadi::DM::reshape(
                    x(Slice(offset, offset + value.numel())), value.rows(),
                    value.columns());
            offset += value.numel();
        }
        return out;
    }

    CasOC::VariablesMX expandVariables(const casadi::MX& x) const {
        CasOC::VariablesMX out;
        using casadi::Slice;
        casadi_int offset = 0;
        for (const auto& key : getSortedVarKeys(m_vars)) {
            const auto& value = m_vars.at(key);
            // Convert a portion of the column vector into a matrix.
            out[key] = casadi::MX::reshape(
                    x(Slice(offset, offset + value.numel())), value.rows(),
                    value.columns());
            offset += value.numel();
        }
        return out;
    }

    /// Flatten the constraints into a row vector, keeping constraints
    /// grouped together by time. Organizing the sparsity of the Jacobian
    /// this way might have benefits for sparse linear algebra.
    template <typename T>
    T flattenConstraints(const Constraints<T>& constraints) const {
        T flat = T(casadi::Sparsity::dense(m_numConstraints, 1));

        int iflat = 0;
        auto copyColumn = [&flat, &iflat](const T& matrix, int columnIndex) {
            using casadi::Slice;
            if (matrix.rows()) {
                flat(Slice(iflat, iflat + matrix.rows())) =
                        matrix(Slice(), columnIndex);
                iflat += matrix.rows();
            }
        };

        // Trapezoidal sparsity pattern (mapping between flattened and expanded
        // constraints) for mesh intervals 0, 1 and 2: Endpoint constraints
        // depend on all time points through their integral.
        //                   0    1    2    3
        //    endpoint       x    x    x    x
        //    path_0         x
        //    kinematic_0    x
        //    residual_0     x
        //    defect_0       x    x
        //    kinematic_1         x
        //    path_1              x
        //    residual_1          x
        //    defect_1            x    x
        //    kinematic_2              x
        //    path_2                   x
        //    residual_2               x
        //    kinematic_3                   x
        //    path_3                        x
        //    residual_3                    x

        // Hermite-Simpson sparsity pattern for mesh intervals 0, 1 and 2:
        //                   0    0.5    1    1.5    2    2.5    3
        //    endpoint       x     x     x     x     x     x     x
        //    path_0         x
        //    kinematic_0    x
        //    residual_0     x
        //    residual_0.5         x
        //    defect_0       x     x     x
        //    interp_con_0   x     x     x
        //    kinematic_1                x
        //    path_1                     x
        //    residual_1                 x
        //    residual_1.5                     x
        //    defect_1                   x     x     x
        //    interp_con_1               x     x     x
        //    kinematic_2                            x
        //    path_2                                 x
        //    residual_2                             x
        //    residual_2.5                                 x
        //    defect_2                               x     x     x
        //    interp_con_2                           x     x     x
        //    kinematic_3                                        x
        //    path_3                                             x
        //    residual_3                                         x
        //                   0    0.5    1    1.5    2    2.5    3

        for (const auto& endpoint : constraints.endpoint) {
            copyColumn(endpoint, 0);
        }

        int igrid = 0;
        // Index for pointsForInterpControls.
        int icon = 0;
        for (int imesh = 0; imesh < m_numMeshPoints; ++imesh) {
            copyColumn(constraints.kinematic, imesh);
            for (const auto& path : constraints.path) {
                copyColumn(path, imesh);
            }
            if (imesh < m_numMeshIntervals) {
                while (m_grid(igrid).scalar() < m_solver.getMesh()[imesh + 1]) {
                    copyColumn(constraints.multibody_residuals, igrid);
                    copyColumn(constraints.auxiliary_residuals, igrid);
                    ++igrid;
                }
                copyColumn(constraints.defects, imesh);
                while (icon < m_pointsForInterpControls.numel() &&
                        m_pointsForInterpControls(icon).scalar() <
                                m_solver.getMesh()[imesh + 1]) {
                    copyColumn(constraints.interp_controls, icon);
                    ++icon;
                }
            }
        }
        // The loop above does not handle the residual at the final grid point.
        copyColumn(constraints.multibody_residuals, m_numGridPoints - 1);
        copyColumn(constraints.auxiliary_residuals, m_numGridPoints - 1);

        OPENSIM_THROW_IF(iflat != m_numConstraints, OpenSim::Exception,
                "Internal error: final value of the index into the flattened "
                "constraints should be equal to the number of constraints.");
        return flat;
    }

    // Expand constraints that have been flattened into a Constraints struct.
    template <typename T>
    Constraints<T> expandConstraints(const T& flat) const {
        using casadi::Sparsity;

        // Allocate memory.
        auto init = [](int numRows, int numColumns) {
            return T(casadi::Sparsity::dense(numRows, numColumns));
        };
        Constraints<T> out;
        out.defects = init(m_numDefectsPerMeshInterval, m_numMeshPoints - 1);
        out.multibody_residuals = init(m_numMultibodyResiduals, 
                m_numGridPoints);
        out.auxiliary_residuals = init(m_numAuxiliaryResiduals,
                m_numGridPoints);
        out.kinematic = init(m_problem.getNumKinematicConstraintEquations(),
                m_numMeshPoints);
        out.endpoint.resize(m_problem.getEndpointConstraintInfos().size());
        for (int iec = 0; iec < (int)m_constraints.endpoint.size(); ++iec) {
            const auto& info = m_problem.getEndpointConstraintInfos()[iec];
            out.endpoint[iec] = init(info.num_outputs, 1);
        }
        out.path.resize(m_problem.getPathConstraintInfos().size());
        for (int ipc = 0; ipc < (int)m_constraints.path.size(); ++ipc) {
            const auto& info = m_problem.getPathConstraintInfos()[ipc];
            out.path[ipc] = init(info.size(), m_numMeshPoints);
        }
        out.interp_controls = init(m_problem.getNumControls(),
                (int)m_pointsForInterpControls.numel());

        int iflat = 0;
        auto copyColumn = [&flat, &iflat](T& matrix, int columnIndex) {
            using casadi::Slice;
            if (matrix.rows()) {
                matrix(Slice(), columnIndex) =
                        flat(Slice(iflat, iflat + matrix.rows()));
                iflat += matrix.rows();
            }
        };

        for (auto& endpoint : out.endpoint) {
            copyColumn(endpoint, 0);
        }

        int igrid = 0;
        // Index for pointsForInterpControls.
        int icon = 0;
        for (int imesh = 0; imesh < m_numMeshPoints; ++imesh) {
            copyColumn(out.kinematic, imesh);
            for (auto& path : out.path) { copyColumn(path, imesh); }
            if (imesh < m_numMeshIntervals) {
                while (m_grid(igrid).scalar() < m_solver.getMesh()[imesh + 1]) {
                    copyColumn(out.multibody_residuals, igrid);
                    copyColumn(out.auxiliary_residuals, igrid);
                    ++igrid;
                }
                copyColumn(out.defects, imesh);
                while (icon < m_pointsForInterpControls.numel() &&
                        m_pointsForInterpControls(icon).scalar() <
                                m_solver.getMesh()[imesh + 1]) {
                    copyColumn(out.interp_controls, icon);
                    ++icon;
                }
            }
        }
        // The loop above does not handle residuals at the final grid point.
        copyColumn(out.multibody_residuals, m_numGridPoints - 1);
        copyColumn(out.auxiliary_residuals, m_numGridPoints - 1);

        OPENSIM_THROW_IF(iflat != m_numConstraints, OpenSim::Exception,
                "Internal error: final value of the index into the flattened "
                "constraints should be equal to the number of constraints.");
        return out;
    }

    ObjectiveBreakdown expandObjectiveTerms(const casadi::DM& terms) const {
        ObjectiveBreakdown out;
        for (int io = 0; io < (int)m_objectiveTermNames.size(); ++io) {
            out.push_back(std::make_pair(
                    m_objectiveTermNames[io], terms(io).scalar()));
        }
        return out;
    }

    friend class NlpsolCallback;
};

} // namespace CasOC

#endif // OPENSIM_CASOCTRANSCRIPTION_H
