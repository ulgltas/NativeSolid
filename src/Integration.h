#pragma once

#include "Config.h"
#include "Structure.h"
#include "Solver.h"

#define PI 3.14159265

class Integration
{
    double totTime;
    double deltaT;
    unsigned long ExtIter;
    std::string algo;
    Solver *solver;

public:
    Integration(Config *config, Structure *structure);
    ~Integration();
    Solver *GetSolver();
    double GettotTime();
    double GetdeltaT();
    void SetExtIter(unsigned long val_ExtIter);
    unsigned long GetExtIter();
    void SetInitialConditions(Config *config, Structure *structure);
    void TemporalIteration(double &t0, double &tf, Structure *structure);
    void StaticIteration(Structure *structure);
    void HarmonicIteration(Config *config, Structure *structure);
    void UpdateSolution();
};
