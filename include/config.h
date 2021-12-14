#pragma once

#include <string>


class Config{

public:
    Config(std::string filename);
    virtual ~Config();
    virtual Config* GetAddress();
    virtual void ReadConfig();
    virtual std::string GetMeshFile();
    virtual std::string GetUnsteady();
    virtual std::string GetCSDSolver();
    virtual std::string GetStructType();
    virtual std::string GetLinearize();
    virtual std::string GetIntegrationAlgo();
    virtual std::string GetRestartSol();
    virtual std::string GetRestartFile();
    virtual std::string GetMovingMarker();
    virtual double GetSpringStiffness();
    virtual double GetSpringMass();
    virtual double GetInertiaCG();
    virtual double GetInertiaFlexural();
    virtual double GetSpringDamping();
    virtual double GetTorsionalStiffness();
    virtual double GetTorsionalDamping();
    virtual double GetCord();
    virtual double GetFlexuralAxis();
    virtual double GetGravityCenter();
    virtual double GetInitialDisp();
    virtual double GetInitialAngle();
    virtual double GetStartTime();
    virtual double GetDeltaT();
    virtual unsigned long GetDeltaIterWrite();
    virtual double GetStopTime();
    virtual double GetOmega();
    virtual unsigned long GetNumberHarmonics();
    virtual double GetRho();

protected:
    std::string ConfigFileName;
    std::string MESH_FILE, UNSTEADY_SIMULATION, CSD_SOLVER, STRUCT_TYPE, LINEARIZE, INTEGRATION_ALGO, RESTART_SOL, RESTART_FILE, MOVING_MARKER;
    double SPRING_STIFFNESS, SPRING_MASS, INERTIA_CG, INERTIA_FLEXURAL, SPRING_DAMPING, TORSIONAL_STIFFNESS, TORSIONAL_DAMPING, CORD, FLEXURAL_AXIS, GRAVITY_CENTER, INITIAL_DISP, INITIAL_ANGLE, START_TIME, DELTA_T, STOP_TIME, OMEGA, RHO;
    unsigned long DELTAITERWRITE, NUMBER_HARMONICS;

};
