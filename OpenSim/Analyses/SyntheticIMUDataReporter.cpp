/* -------------------------------------------------------------------------- *
 *                      OpenSim:  SyntheticIMUDataReporter.cpp                *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2020 Stanford University and the Authors                *
 * Author(s): Ayman Habib                                                     *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */


//=============================================================================
// INCLUDES
//=============================================================================
#include <OpenSim/Common/IO.h>
#include <OpenSim/Common/STOFileAdapter.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Analyses/SyntheticIMUDataReporter.h>
#include <OpenSim/Simulation/OpenSense/SyntheticIMU.h>

using namespace OpenSim;
using namespace std;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Destructor.
 */
SyntheticIMUDataReporter::~SyntheticIMUDataReporter()
{
}
//_____________________________________________________________________________
/**
 */
SyntheticIMUDataReporter::SyntheticIMUDataReporter(Model *aModel) :
    Analysis(aModel), _modelLocal(nullptr) {
    setNull();
    constructProperties();
    if(aModel) 
        setModel(*aModel);
}
// Copy constructor and virtual copy 
//_____________________________________________________________________________
/**
 * Copy constructor.
 *
 */
SyntheticIMUDataReporter::SyntheticIMUDataReporter(const SyntheticIMUDataReporter &aSyntheticIMUDataReporter): 
    Analysis(aSyntheticIMUDataReporter),
    _modelLocal(nullptr) {
    setNull();
    // COPY TYPE AND NAME
    *this = aSyntheticIMUDataReporter;
}

//=============================================================================
// OPERATORS
//=============================================================================
//-----------------------------------------------------------------------------
// ASSIGNMENT
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
/**
 * Assign this object to the values of another.
 *
 * @return Reference to this object.
 */
SyntheticIMUDataReporter& SyntheticIMUDataReporter::
operator=(const SyntheticIMUDataReporter& other)
{
    // BASE CLASS
    Analysis::operator=(other);
    copyProperty_report_orientations(other);
    copyProperty_report_angular_velocities(other);
    copyProperty_report_linear_accelerations(other);
    copyProperty_IMU_frames(other);
    copyProperty_frame_paths(other);
    _modelLocal = nullptr;
    return(*this);
}

//_____________________________________________________________________________
/**
 * SetNull().
 */
void SyntheticIMUDataReporter::setNull()
{
    setAuthors("Ayman Habib");

    setName("SyntheticIMUDataReporter");
    _modelLocal = nullptr;
}

//=============================================================================
// CONSTRUCTION METHODS
//=============================================================================
//=============================================================================
// GET AND SET
//=============================================================================
//_____________________________________________________________________________
//=============================================================================
// ANALYSIS
//=============================================================================
//_____________________________________________________________________________
/**
 * Record the results.
 */
int SyntheticIMUDataReporter::
record(const SimTK::State& s)
{
    if(_modelLocal == nullptr) return -1;

    // Set model to whatever defaults have been updated to from the last iteration
    SimTK::State& sWorkingCopy = _modelLocal->updWorkingState();
    sWorkingCopy.setTime(s.getTime());
    //_modelLocal->initStateWithoutRecreatingSystem(sWorkingCopy); 

    // update Q's and U's
    sWorkingCopy.setQ(s.getQ());
    sWorkingCopy.setU(s.getU());

    _modelLocal->getMultibodySystem().realize(
            sWorkingCopy, SimTK::Stage::Report);

    return 0;
}
//_____________________________________________________________________________
/**
 * This method is called at the beginning of an analysis so that any
 * necessary initializations may be performed.
 *
 * This method is meant to be called at the beginning of an integration 
 *
 * @param s Current state .
 *
 * @return -1 on error, 0 otherwise.
 */
int SyntheticIMUDataReporter::begin(const SimTK::State& s )
{
    if(!proceed()) return(0);

    _orientationsReporter.clearTable();
    _angularVelocityReporter.clearTable();
    _linearAccelerationsReporter.clearTable();

    if (_imuComponents.empty()) {
        // Populate _imuComponents based on properties
        if (get_IMU_frames() == "Bodies") { 
            reportAllBodies();
        } else if (get_IMU_frames() == "Frames") {
            reportAllFrames();
        } else if (get_IMU_frames() == "Custom") {
            for (int i = 0; i < getProperty_frame_paths().size(); ++i) {
                // This may need stricter error checking to make sure we have a valid Frame/Path
                _imuComponents.push_back(ComponentPath(get_frame_paths(i)));
            }
        } else {
            log_warn("SyntheticIMUDataReporter has invalid specification");
            log_warn("Current selection is {}, required Bodies/Frames/Custom", get_IMU_frames());
            log_warn("All bodies will be reported.");
            reportAllBodies();
        }
    }
    // If alreay part of the system, then a rerun and no need to add to _modelLocal
    if (!_orientationsReporter.hasSystem()) {
        _modelLocal->addComponent(&_orientationsReporter);
        _modelLocal->addComponent(&_angularVelocityReporter);
        _modelLocal->addComponent(&_linearAccelerationsReporter);

        for (auto& path : _imuComponents) {
            const Component& comp = _modelLocal->getComponent(path.toString());
            if (get_report_orientations())
                _orientationsReporter.addToReport(
                    comp.getOutput("rotation_as_quaternion"), comp.getName());
            if (get_report_angular_velocities())
                _angularVelocityReporter.addToReport(
                    comp.getOutput("angular_velocity"), comp.getName());
            if (get_report_linear_accelerations())
                _linearAccelerationsReporter.addToReport(
                    comp.getOutput("linear_acceleration"), comp.getName());
        }
    }
    _modelLocal->initSystem();

    // RECORD
    int status = record(s);
 
    return(status);
}
//_____________________________________________________________________________
/**
 * This method is called to perform the analysis.  It can be called during
 * the execution of a forward integrations or after the integration by
 * feeding it the necessary data.
 *
 * This method should be overridden in derived classes.  It is
 * included here so that the derived class will not have to implement it if
 * it is not necessary.
 *
 * @param s Current state .
 *
 * @return -1 on error, 0 otherwise.
 */
int SyntheticIMUDataReporter::step(const SimTK::State& s, int stepNumber )
{
    if(!proceed(stepNumber)) return(0);

    record(s);

    return(0);
}
//_____________________________________________________________________________
/**
 * This method is called at the end of an analysis so that any
 * necessary finalizations may be performed.
 *
 * @param s Current state 
 *
 * @return -1 on error, 0 otherwise.
 */
int SyntheticIMUDataReporter::end( const SimTK::State& s )
{
    if(!proceed()) return(0);

    record(s);

    return(0);
}


//=============================================================================
// IO
//=============================================================================
//_____________________________________________________________________________
/**
 * Print results.
 * 
 * The file names are constructed as
 * aDir + "/" + aBaseName + "_" + ComponentName + aExtension
 *
 * @param aDir Directory in which the results reside.
 * @param aBaseName Base file name.
 * @param aDT Desired time interval between adjacent storage vectors.  Linear
 * interpolation is used to print the data out at the desired interval.
 * @param aExtension File extension.
 *
 * @return 0 on success, -1 on error.
 */
int SyntheticIMUDataReporter::
printResults(const string &aBaseName,const string &aDir,double aDT,
                 const string &aExtension)
{
    auto& rotationsTable = _orientationsReporter.getTable();
    auto& angVelTable = _angularVelocityReporter.getTable();
    auto& linAccTable = _linearAccelerationsReporter.getTable();
    {
        IO::CwdChanger cwd = IO::CwdChanger::changeTo(aDir);
        if (get_report_orientations())
            STOFileAdapter_<SimTK::Quaternion>::write(
                rotationsTable, aBaseName +"_"+"orientations.sto");
        if (get_report_angular_velocities())
            STOFileAdapter_<SimTK::Vec3>::write(
                angVelTable, aBaseName + "_" + "angular_veolcity.sto");
        if (get_report_linear_accelerations())
            STOFileAdapter_<SimTK::Vec3>::write(
                linAccTable, aBaseName + "_" + "linear_accelerations.sto");
    }
    return(0);
}
template <typename T>
void SyntheticIMUDataReporter::reportAll() {
    _modelLocal = std::make_unique<Model>();
    _modelLocal.reset(_model->clone());
    ComponentList <T> frameComponents =
            _modelLocal->getComponentList<T>();
    // Add SyntheticIMUs to the model if needed
    for (auto& frame : frameComponents) {
        // Ceate Synthetic IMU, connect it to frame
        SyntheticIMU* next_imu = new SyntheticIMU();
        next_imu->setName(frame.getName() + "_imu");
        next_imu->connectSocket_attachment_frame(frame);
        _modelLocal->addComponent(next_imu);
        _imuComponents.push_back(next_imu->getAbsolutePath());
    };
    _modelLocal->finalizeConnections();
}