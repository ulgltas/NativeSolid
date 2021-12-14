#include "../include/integration.h"
#include "../include/solver.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <cmath>

using namespace std;

Integration::Integration(Config *config, Structure *structure){

  solver = NULL;
  bool linear;

  linear = (config->GetLinearize()) == "YES";

  if (config->GetUnsteady() == "YES") {
    if(config->GetIntegrationAlgo() == "ALPHAGEN"){
      solver = new AlphaGenSolver(structure->GetnDof(), config->GetRho(), linear);
    }
    else if(config->GetIntegrationAlgo() == "RK4"){
      solver = new RK4Solver(structure->GetnDof(), linear);
    }
    else{
    }
  }
  else if (config->GetUnsteady() == "HARMONIC"){
    solver = new HarmonicSolver(structure->GetnDof(), config->GetNumberHarmonics(), linear);
  }
  else
    solver = new StaticSolver(structure->GetnDof(), linear);
}

Integration::~Integration(){

  if (solver != NULL) delete solver;
}

Solver* Integration::GetSolver(){

  return solver;
}

double Integration::GettotTime(){
  return totTime;
}

double Integration::GetdeltaT(){
  return deltaT;
}

void Integration::SetExtIter(unsigned long val_ExtIter){

    ExtIter = val_ExtIter;
}

unsigned long Integration::GetExtIter(){

    return ExtIter;
}

void Integration::SetInitialConditions(Config* config, Structure* structure){

  ExtIter = 0;
  solver->SetInitialState(config, structure);
  structure->SetCenterOfRotation_Y((solver->GetDisp())[0]);

}

void Integration::TemporalIteration(double& t0, double& tf, Structure *structure){

  int rank = 0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  double currentTime(tf);
  deltaT = tf - t0;
  solver->Iterate(t0, tf, structure);

  if(rank == 0){
  if(structure->GetnDof() == 1){
    cout << fixed
         << setw(10) << "Time"
         << setw(15) << "Displacement"
         << setw(15) << "Velocity"
         << setw(15) << "Acceleration" << endl;
    cout << fixed
         << setw(10) << currentTime
         << setw(15) << (solver->GetDisp())[0]
         << setw(15) << (solver->GetVel())[0]
         << setw(15) << (solver->GetAcc())[0] << endl;
  }
  else if(structure->GetnDof() == 2){cout << fixed
         << setw(10) << "Time"
         << setw(15) << "Displacement1"
         << setw(15) << "Displacement2"
         << setw(15) << "Velocity1"
         << setw(15) << "Velocity2"
         << setw(15) << "Acceleration1"
         << setw(15) << "Acceleration2" << endl;
    cout << fixed
         << setw(10) << currentTime
         << setw(15) << (solver->GetDisp())[0]
         << setw(15) << (solver->GetDisp())[1]
         << setw(15) << (solver->GetVel())[0]
         << setw(15) << (solver->GetVel())[1]
         << setw(15) << (solver->GetAcc())[0]
         << setw(15) << (solver->GetAcc())[1] << endl;}
  }
}


void Integration::StaticIteration(Structure *structure){

  int rank = 0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  double t0 = 0.;
  double tf = 0.;
  deltaT = 0.;
  solver->Iterate(t0, tf, structure);

  if(rank == 0){
    if(structure->GetnDof() == 1){
      cout << fixed
           << setw(10) << "Time"
           << setw(15) << "Displacement" << endl;
      cout << fixed
           << setw(10) << tf
           << setw(15) << (solver->GetDisp())[0] << endl;
    }
    else if(structure->GetnDof() == 2){
      cout << fixed
           << setw(10) << "Time"
           << setw(15) << "Displacement1"
           << setw(15) << "Displacement2" << endl;
      cout << fixed
           << setw(10) << tf
           << setw(15) << (solver->GetDisp())[0]
           << setw(15) << (solver->GetDisp())[1] << endl;
    }
  }
}

void Integration::HarmonicIteration(Config* config, Structure *structure){

  int rank = 0;
  double currentTime = 0.0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  double t0 = 0.;
  double tf = 0.;
  unsigned long nInstances;
  deltaT = config->GetDeltaT();
  solver->Iterate(t0, deltaT, structure);
  nInstances = 2*config->GetNumberHarmonics()+1;
  
  if(rank == 0){
  if(structure->GetnDof() == 1){
    cout << fixed
         << setw(10) << "Time"
         << setw(15) << "Displacement"
         << setw(15) << "Velocity"
         << setw(15) << "Acceleration" << endl;
    for (int i = 0; i < nInstances; i++)
    {
      cout << fixed
         << setw(10) << currentTime
         << setw(15) << (solver->GetDisp())[i]
         << setw(15) << (solver->GetVel())[i]
         << setw(15) << (solver->GetAcc())[i] << endl;
      currentTime += config->GetStopTime()/nInstances;
    }
    
    
  }
  else if(structure->GetnDof() == 2){
    cout << fixed
         << setw(10) << "Time"
         << setw(15) << "Displacement1"
         << setw(15) << "Displacement2"
         << setw(15) << "Velocity1"
         << setw(15) << "Velocity2"
         << setw(15) << "Acceleration1"
         << setw(15) << "Acceleration2" << endl;
    currentTime = 0.;
    for (int i = 0; i < nInstances; i++){
        cout << fixed
              << setw(10) << currentTime
              << setw(15) << (solver->GetDisp())[i]
              << setw(15) << (solver->GetDisp())[i+nInstances]
              << setw(15) << (solver->GetVel())[i]
              << setw(15) << (solver->GetVel())[i+nInstances]
              << setw(15) << (solver->GetAcc())[i]
              << setw(15) << (solver->GetAcc())[i+nInstances] << endl;
        currentTime += config->GetStopTime()/nInstances;
      }
    }
  }
}

void Integration::UpdateSolution(){

  solver->SaveToThePast();

}
