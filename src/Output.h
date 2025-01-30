#pragma once

#include <fstream>

#include "MatVec.h"
#include "Config.h"
#include "Structure.h"
#include "Integration.h"

class Output
{
public:
    Output();
    // void WriteHistory(Integration *solver, Structure *structure,
    //                   std::ofstream *outputfile, const double &time);
    // void WriteRestart(Integration *solver, Structure *structure);
    void WriteRestart(Integration *solver, Structure *structure, Config *config);
    void WriteAdjointOutput(double dampingDer);
    // void WriteStaticSolution(Config *config, Integration *solver,
    //                          Structure *structure, std::ofstream *outputfile);
    // void WriteRestart();
};
