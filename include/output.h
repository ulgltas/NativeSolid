#pragma once

#include <fstream>

#include "MatVec.h"
#include "config.h"
#include "structure.h"
#include "integration.h"

class Output{

public:
    Output(void);
    ~Output();
    //void WriteHistory(Integration* solver, Structure* structure, std::ofstream* outputfile, const double & time);
    void WriteRestart(Integration* solver, Structure* structure, Config* config);
    //void WriteStaticSolution(Config* config, Integration* solver, Structure* structure, std::ofstream* outputfile);
    //void WriteRestart();
};
