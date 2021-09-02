#ifndef _OPENSIM_DC_MOTOR_H_
#define _OPENSIM_DC_MOTOR_H_
/* -------------------------------------------------------------------------- *
*                             DCMotor.h                                *
* -------------------------------------------------------------------------- *
* The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
* See http://opensim.stanford.edu and the NOTICE file for more information.  *
* OpenSim is developed at Stanford University and supported by the US        *
* National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
* through the Warrior Web program.                                           *
*                                                                            *
* Copyright (c) 2005-2012 Stanford University                                *
*                                                                            *                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License"); you may    *
* not use this file except in compliance with the License. You may obtain a  *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied    *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
* -------------------------------------------------------------------------- */

/**
* A class that implements an DC motor that applies equal and opposite torques
* on the two bodies: rotor and stator.
* The model of the motor contain only the DC torque equation. The mass and
* inertia of rotor must be included as part of system's mass. The torque is applied about
* the shaft axis specified in global by default, otherwise in rotor's coordinate system.
*
* Author: Vinh Q. Nguyen, Andrew Lapre, Mark Price - Umass Amherst
* Contact: vinhnguyen@umass.edu
* Version 1.0 (Nov 2016)
* This class is based on the template of 'TorqueActuator' class.
*/

#include <OpenSim/Actuators/osimActuatorsDLL.h>

#include <string>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Force.h>
#include <OpenSim/Simulation/Model/Actuator.h>

namespace OpenSim {

	class PhysicalFrame;
	class Model;

	class OSIMACTUATORS_API DCMotor : public ScalarActuator {
		OpenSim_DECLARE_CONCRETE_OBJECT(DCMotor, ScalarActuator);


		//==============================================================================================
		//											 PUBLIC
		//==============================================================================================
	public:
		//--------------------------------------------------------------------------------
		// Declare parameters and parts of motor
		//--------------------------------------------------------------------------------

		OpenSim_DECLARE_OPTIONAL_PROPERTY(rotor, std::string,
			"Name of Body Rotor to which the torque of motor is applied.");
		OpenSim_DECLARE_OPTIONAL_PROPERTY(stator, std::string,
			"Name of Body Stator to which the equal and opposite torque is applied.");
		/* The default is torque_is_global=true. */
		OpenSim_DECLARE_PROPERTY(axis_is_global, bool,
			"Interpret axis in Ground frame if true; otherwise, Rotor's frame.");
		/* The default direction for the axis is z (0,0,1). */
		OpenSim_DECLARE_PROPERTY(shaft_axis, SimTK::Vec3,
			"Fixed direction about which torque is applied, in Ground or Rotor "
                "frame depending on 'torque_is_global' property.");
        /* The default is force_current_limit=false. */
        OpenSim_DECLARE_PROPERTY(enforce_current_limit, bool,
                "Enforce current limit in dynamic equations (bad for direct collocation).");

		OpenSim_DECLARE_PROPERTY(optimal_force, double,
			"The maximum torque produced by this actuator when fully activated.");

		/* Parameters for Motor. */
		OpenSim_DECLARE_PROPERTY(K, double, "Motor torque constant (Nm/A) (default value: 0.0)");

		OpenSim_DECLARE_PROPERTY(L, double, "Motor inductance (H) (default value: 0.0)");

		OpenSim_DECLARE_PROPERTY(R, double, "Motor resistance (Ohm) (default value: 0.0)");

		OpenSim_DECLARE_PROPERTY(b, double, "Rotor Damping (default value: 0.0)"
			"If specifing here, then there should be no more damping when defining the joint between rotor and stator");

		OpenSim_DECLARE_PROPERTY(maximum_current, double, "maximum_current (A) (default value: 10)");

		OpenSim_DECLARE_PROPERTY(initial_current, double, "initial_current (A) (default value: 0.0)");

		OpenSim_DECLARE_PROPERTY(voltageRatio, double, "voltage control scale factor (default value: 1.0)");


		/* Default constructor leaves body names unspecified. */
		DCMotor();
		/* Destructor. */
		~DCMotor();


		//--------------------------------------------------------------------------
		// Get and Set Motor Parameters and Properties
		//--------------------------------------------------------------------------
		/* Set the 'shaft_axis' property; frame is interpreted
		according to the 'axis_is_global' property. **/
		void setShaftAxis(const SimTK::Vec3& axis)
		{
			set_shaft_axis(axis);
		}
		/* Return the current value of the 'shaft_axis' property. */
		const SimTK::Vec3& getShaftAxis() const
		{
			return get_shaft_axis();
		}

		/* Set the 'axis_is_global' property that determines how to interpret
		the 'shaft_axis' vector; if not global (Ground frame) it is in rotor's frame. */
		void setAxisIsGlobal(bool isGlobal)
		{
			set_axis_is_global(isGlobal);
		}
		/* Return the current value of the 'torque_is_global' property. */
		bool getAxisIsGlobal() const
		{
			return get_axis_is_global();
		}

		void setVoltageRatio(double volRatio)
		{
			set_voltageRatio(volRatio);
		}

		double getVoltageRatio() const
		{
			return get_voltageRatio();
		}

		/** Set the 'optimal_force' property. **/
		void setOptimalForce(double optimalForce)
		{
			set_optimal_force(optimalForce);
		}
		/** Get the current value of the 'optimal_force' property. **/
		double getOptimalForce() const // Part of Actuator interface.
		{
			return get_optimal_force();
		}


		/* Motor parameters. */
		double getK() const;  void setK(double K);

		double getL() const;  void setL(double L);

		double getR() const;  void setR(double R);

		double getRotorDamping() const; void setRotorDamping(double b);

		/* Limit the power supplied. */
		double getMaximumCurrent() const; void setMaximumCurrent(double maximum_current);

		/* Set the bodies (stator and rotor) . */
		void setRotor(const PhysicalFrame& body);
		void setStator(const PhysicalFrame& body);

		/* Get the bodies (stator and rotor). */
		const PhysicalFrame& getRotor() const { return *_rotor; }
		const PhysicalFrame& getStator() const { return *_stator; }


		//==============================================================================================
		//									    PRIVATE
		//==============================================================================================
	private:
		void constructProperties();

		//--------------------------------------------------------------------------
		// Implement Force interface
		//--------------------------------------------------------------------------
		void computeForce(const SimTK::State& state,
			SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
			SimTK::Vector& mobilityForces) const override;

		//--------------------------------------------------------------------------
		// Implement Actuator interface
		//--------------------------------------------------------------------------
		double computeActuation(const SimTK::State& s) const override;
		// Return the stress, defined as abs(torque/Maximum torque).
		double getStress(const SimTK::State& state) const override;

		//--------------------------------------------------------------------------
		// Implement ModelComponent interface
		//--------------------------------------------------------------------------
		void extendAddToSystem(SimTK::MultibodySystem& system) const override;
		void extendInitStateFromProperties(SimTK::State& state) const override;
		void extendSetPropertiesFromState(const SimTK::State& state) override;

		double getCurrent(const SimTK::State& s) const;
		void setCurrent(SimTK::State& s, double current) const;

		double getCurrentDeriv(const SimTK::State& s) const;
		void setCurrentDeriv(const SimTK::State& s, double activeMotorUnitsDeriv) const;

		//void setStateVariableDeriv(const SimTK::State& s, const std::string &aStateName, double aValue) const;
		//double getStateVariableDeriv(const SimTK::State& s, const std::string &aStateName) const;

		void computeStateVariableDerivatives(const SimTK::State& s) const override;

		// Setup method initializes Body reference pointers to match the names.
		void extendConnectToModel(Model& model) override;

		//--------------------------------------------------------------------------
		// Implement Object interface.
		//--------------------------------------------------------------------------
		void updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber = -1);

		//--------------------------------------------------------------------------
		// Data
		//--------------------------------------------------------------------------
		// Note: reference pointers are automatically set to null on construction 
		// and also on copy construction and copy assignment.

		// Corresponding Body to which the torque actuator is applied.
		SimTK::ReferencePtr<const PhysicalFrame> _rotor;

		// Corresponding Body to which the equal and opposite torque is applied.
		SimTK::ReferencePtr<const PhysicalFrame> _stator;

		//==============================================================================
	};	// END of class DCMotor


} // end of namespace OpenSim

#endif // _OPENSIM_DC_MOTOR_H_


