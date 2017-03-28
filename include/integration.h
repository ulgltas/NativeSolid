#pragma once

//#include "MatVec.h"
#include "config.h"
#include "structure.h"
#include "solver.h"

#define PI 3.14159265


class Integration{

protected:
    double totTime;
    double deltaT;
    unsigned long ExtIter;
    std::string algo;

    Solver* solver;

public:
    Integration(Config* config, Structure* structure);
    ~Integration();
    Solver* GetSolver();
    double GettotTime();
    double GetdeltaT();
    void SetExtIter(unsigned long val_ExtIter);
    unsigned long GetExtIter();
    //void SetIntegrationParam(Config* config);
    //void SetLoadsAtTime(Config* config, Structure* structure, const double & time, double FSI_Load);
    void SetInitialConditions(Config* config, Structure* structure);
    void TemporalIteration(double& t0, double& tf, Structure *structure);
    //void StaticIteration(Config *config, Structure *structure);
    void UpdateSolution();
};
