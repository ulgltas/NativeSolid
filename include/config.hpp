#pragma once

#ifdef HAVE_MPI
    #include "mpi.h"
#endif

#include <cstring>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <algorithm>

using namespace std;

class Config{

public:
    Config(string filename);
    virtual ~Config();
    virtual Config* GetAddress();
    virtual void ReadConfig();
    virtual string GetUnsteady();
    virtual string GetCSDSolver();
    virtual string GetStructType();
    virtual string GetForceInputType();
    virtual string GetForceInputFile();
    virtual string GetAnalyticalFunction();
    virtual string GetIntegrationAlgo();
    virtual string GetRestartSol();
    virtual string GetRestartFile();
    virtual double GetSpringStiffness();
    virtual double GetSpringMass();
    virtual double GetSpringDamping();
    virtual double GetTorsionalStiffness();
    virtual double GetTorsionalDamping();
    virtual double GetCord();
    virtual double GetFlexuralAxis();
    virtual double GetStartTime();
    virtual double GetDeltaT();
    virtual double GetStopTime();
    virtual double GetFrequency();
    virtual double GetAmplitude();
    virtual double GetConstantForceValue();
    virtual double GetAlpha();
    virtual double GetRho();

protected:
    string ConfigFileName;
    string UNSTEADY_SIMULATION, CSD_SOLVER, STRUCT_TYPE, FORCE_INPUT_TYPE, FORCE_INPUT_FILE, ANALYTICAL_FUNCTION, INTEGRATION_ALGO, RESTART_SOL, RESTART_FILE;
    double SPRING_STIFFNESS, SPRING_MASS, SPRING_DAMPING, TORSIONAL_STIFFNESS, TORSIONAL_DAMPING, CORD, FLEXURAL_AXIS, START_TIME, DELTA_T, STOP_TIME, FREQUENCY, AMPLITUDE, CONSTANT_FORCE_VALUE, ALPHA, RHO;

};
