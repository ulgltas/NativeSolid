#include "Output.h"
#include <iostream>

Output::Output() {}

/*
void Output::WriteHistory(Integration* solver, Structure* structure, std::ofstream* outputfile, const double & time){

  if(structure->GetnDof() == 1){
    if(time == 0){
      cout << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << std::endl;
      outputfile[0] << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << std::endl;
    }
    cout << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << std::endl;
    outputfile[0] << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << std::endl;
  }
  else if(structure->GetnDof() == 2){
    if(time == 0){
      cout << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << std::endl;
      outputfile[0] << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << std::endl;
    }
    cout << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << std::endl;
    outputfile[0] << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << std::endl;
  }
}*/

void Output::WriteRestart(Integration *integrator, Structure *structure, Config *config)
{
    std::ofstream RestartFile;
    double time(0.0);
    unsigned short nInst = 2 * config->GetNumberHarmonics() + 1;
    RestartFile.open("restart_solid.dat", std::ios::out);
    RestartFile.precision(8);

    if (structure->GetnDof() == 1)
    {
        RestartFile << "Displacement" << "\t" << "Velocity" << "\t" << "Acceleration" << std::endl;
        for (unsigned short iInst = 0; iInst < nInst; iInst++)
        {
            RestartFile << std::fixed
                        << (integrator->GetSolver()->GetDisp())[iInst] << "\t"
                        << (integrator->GetSolver()->GetVel())[iInst] << "\t"
                        << (integrator->GetSolver()->GetAcc())[iInst] << std::endl;
        }
    }
    else if (structure->GetnDof() == 2)
    {
        RestartFile << "Time" << "\t" << "Displacement_1" << "\t" << "Displacement_2" << "\t" << "Velocity_1" << "\t" << "Velocity_2" << "\t" << "Acceleration_1" << "\t" << "Acceleration_2" << std::endl;
        for (unsigned short iInst = 0; iInst < nInst; iInst++)
        {
            RestartFile << std::fixed
                        << time << "\t"
                        << (integrator->GetSolver()->GetDisp())[iInst] << "\t"
                        << (integrator->GetSolver()->GetDisp())[iInst + nInst] << "\t"
                        << (integrator->GetSolver()->GetVel())[iInst] << "\t"
                        << (integrator->GetSolver()->GetVel())[iInst + nInst] << "\t"
                        << (integrator->GetSolver()->GetAcc())[iInst] << "\t"
                        << (integrator->GetSolver()->GetAcc())[iInst + nInst] << std::endl;
            time += 0.1;
        }
    }
    RestartFile.close();
}

void Output::WriteAdjointOutput(double dampingDer)
{
    std::ofstream AdjointFile;

    AdjointFile.open("solid_gradients.csv", std::ios::out);
    AdjointFile.precision(12);

    AdjointFile << "plunge_damping" << std::endl;
    AdjointFile << dampingDer << std::endl;

    AdjointFile.close();
}

/*void Output::WriteStaticSolution(Config* config, Integration* solver, Structure* structure, std::ofstream* outputfile){
  if(structure->GetnDof() == 1){
    cout << "Static displacement is : " << (*(solver->GetDisp()))[0] << " [m]" << std::endl;
    cout << "Writing displacement into a solution file" << std::endl;
    if(config->GetStructType() == "SPRING_HOR"){
      outputfile[0] << -1.000 << "\t" << (*(solver->GetDisp()))[0] << "\t" << 0.000;
    }
    else if(config->GetStructType() == "SPRING_VER"){
      outputfile[0] << -1.000 << "\t" << 0.000 << "\t" << (*(solver->GetDisp()))[0];
    }
  }
  else if(structure->GetnDof() == 2){
    cout << "Plunging displacement is :" << (*(solver->GetDisp()))[0] << " [m]" << std::endl;
    cout << "Pitching rotation is :" << (*(solver->GetDisp()))[1] << " [rad]" << std::endl;
    cout << "Writing displacement into a solution file" << std::endl;
    outputfile[0] << -1.000 << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1];
  }
}
*/
