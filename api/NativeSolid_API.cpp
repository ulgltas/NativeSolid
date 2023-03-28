#include "NativeSolid_API.h"
#include <cmath>
#include <iomanip>

using namespace std;

const int MASTER_NODE = 0;


NativeSolidSolver::NativeSolidSolver(string str, bool FSIComp):confFile(str){

  int rank = MASTER_NODE;

/*--- MPI initialization, and buffer setting ---*/
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if(rank == MASTER_NODE){
    if(FSIComp)cout << endl << "***************************** Setting NativeSolid for FSI simulation *****************************" << endl;
    else cout << endl <<"***************************** Setting NativeSolid for CSD simulation *****************************" << endl;
  }

  config = NULL;
  output = NULL;
  geometry = NULL;
  structure = NULL;
  integrator = NULL;

  /*--- Initialize the main containers ---*/
  config = new Config(confFile);
  output = new Output();

  /*--- Read CSD configuration file ---*/
  cout << endl << "\n----------------------- Reading Native configuration file ----------------------" << endl;
  config->ReadConfig();

  /*--- Read a SU2 native mesh file ---*/
  cout << endl << "\n----------------------- Reading SU2 based mesh file ----------------------" << endl;
  geometry = new Geometry(config);

  /*--- Initialize structural container and create the structural model ---*/
  cout << endl << "\n----------------------- Creating the structural model ----------------------" << endl;
  structure = new Structure(config);

  cout << endl << "\n----------------------- Creating the FSI interface ----------------------" << endl;
  double* Coord;
  unsigned long iMarker(0), iPoint;

  if (config->GetKindProblem() == "ADJOINT"){
    historyFile.open("NativeHistoryAdjoint.dat", ios::out);
    historyFile2.open("NativeHistoryAdjointFSI.dat", ios::out);
  }
  else{
    historyFile.open("NativeHistory.dat", ios::out);
    historyFile2.open("NativeHistoryFSI.dat", ios::out);
  }
  historyFile.precision(8); // More digits than default for history files
  historyFile2.precision(8);
  

  while(iMarker < geometry->GetnMarkers()){
    cout << iMarker << " " << geometry->GetMarkersMoving(iMarker) << " " << geometry->GetnMarkers() << endl;
    if(geometry->GetMarkersMoving(iMarker)){
      nSolidInterfaceVertex = geometry->nVertex[iMarker];
    break;
    }
    iMarker++;
  }

  varCoordNorm = 0.0;

  cout << nSolidInterfaceVertex << " nodes on the moving interface have to be tracked." << endl;

  /*--- Initialize the temporal integrator ---*/
  cout << endl << "\n----------------------- Setting integration parameter ----------------------" << endl;
  integrator = new Integration(config, structure);
  integrator->SetInitialConditions(config, structure);

  cout << endl << "\n----------------------- Setting FSI features ----------------------" << endl;
  q_uM1.Initialize(structure->GetnDof());
  q_uM1.Reset();

  posDV.Initialize(config->GetNumberDesignVariables());
  posDV.Reset();
  magDV.Initialize(config->GetNumberDesignVariables());
  magDV.Reset();
  sideDV.Initialize(config->GetNumberDesignVariables());
  sideDV.Reset();

  if(rank == MASTER_NODE){
    if(structure->GetnDof() == 1){
      if((config->GetUnsteady() == "YES" || config->GetUnsteady() == "HARMONIC") && config->GetKindProblem() != "ADJOINT"){
        historyFile << fixed
                    << setw(10) << "Delta_t"
                    << setw(10) << "FSI iter"
                    << setw(15) << "Displacement"
                    << setw(15) << "Velocity"
                    << setw(15) << "Acceleration" << endl;
        historyFile2 << fixed
                      << setw(10) << "Time"
                      << setw(10) << "FSI iter"
                      << setw(15) << "Displacement"
                      << setw(15) << "Velocity"
                      << setw(15) << "Acceleration" << endl;
      }
      else if (config->GetKindProblem() == "ADJOINT"){
        historyFile << fixed
                    << setw(10) << "Delta_t"
                    << setw(10) << "FSI iter"
                    << setw(15) << "dE/dkh" << endl;
        historyFile2 << fixed
                      << setw(10) << "FSI iter"
                      << setw(15) << "dE/dkh" << endl;
      }
      else{
        historyFile << fixed
                    << setw(10) << "Delta_t"
                    << setw(10) << "FSI iter"
                    << setw(15) << "Displacement" << endl;
        historyFile2 << fixed
                      << setw(10) << "FSI iter"
                      << setw(15) << "Displacement" << endl;
      }
    }
    else if(structure->GetnDof() == 2){
      if(config->GetUnsteady() == "YES" || config->GetUnsteady() == "HARMONIC" && config->GetKindProblem() != "ADJOINT"){
        historyFile << fixed
                    << setw(10) << "Delta_t"
                    << setw(10) << "FSI iter"
                    << setw(15) << "Displacement 1"
                    << setw(15) << "Displacement 2"
                    << setw(15) << "Velocity 1"
                    << setw(15) << "Velocity 2"
                    << setw(15) << "Acceleration 1"
                    << setw(15) << "Acceleration 2" << endl;
        historyFile2 << fixed
                      << setw(10) << "Time"
                      << setw(10) << "FSI iter"
                      << setw(15) << "Displacement 1"
                      << setw(15) << "Displacement 2"
                      << setw(15) << "Velocity 1"
                      << setw(15) << "Velocity 2"
                      << setw(15) << "Acceleration 1"
                      << setw(15) << "Acceleration 2" << endl;
      }
      else if (config->GetKindProblem() == "ADJOINT"){
        historyFile << fixed
                    << setw(10) << "Delta_t"
                    << setw(10) << "FSI iter"
                    << setw(15) << "Energy"
                    << setw(15) << "dE/dkh"
                    << setw(15) << "dE/dka" << endl;
        historyFile2 << fixed
                      << setw(10) << "FSI iter"
                      << setw(15) << "Energy"
                      << setw(15) << "dE/dkh"
                      << setw(15) << "dE/dka" << endl;
      }
      else{
        historyFile << fixed
                    << setw(10) << "Delta_t"
                    << setw(10) << "FSI iter"
                    << setw(15) << "Displacement 1"
                    << setw(15) << "Displacement 2" << endl;
        historyFile2  << fixed
                      << setw(10) << "FSI iter"
                      << setw(15) << "Displacement 1"
                      << setw(15) << "Displacement 2" << endl;
      }
    }
  }

  if(rank == MASTER_NODE){
    if(FSIComp)cout << endl << "***************************** NativeSolid is set for FSI simulation *****************************" << endl;
    else cout << endl <<"***************************** NativeSolid is set for CSD simulation *****************************" << endl;
  }

}

NativeSolidSolver::~NativeSolidSolver(){}

void NativeSolidSolver::exit(){

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  historyFile.close();
  if(rank == MASTER_NODE){
    cout << "Solid history file is closed." << endl;
    cout << endl << "***************************** Exit NativeSolid *****************************" << endl;
  }

  if (config != NULL) delete config;
  if (geometry != NULL) delete geometry;
  if (structure != NULL) delete structure;
  if (integrator != NULL) delete integrator;
  if (output != NULL) delete output;

}

void NativeSolidSolver::preprocessIteration(unsigned long ExtIter){

    integrator->SetExtIter(ExtIter);
}

void NativeSolidSolver::timeIteration(double t0, double tf){

  integrator->TemporalIteration(t0, tf, structure);

  //mapRigidBodyMotion(false, false);
  computeInterfacePosVel(false);

}

/*
void NativeSolidSolver::mapRigidBodyMotion(bool prediction, bool initialize){

  double *Coord, *Coord_n, newCoord[3], Center[3], Center_n[3], newCenter[3], rotCoord[3], r[3];
  double varCoord[3] = {0.0, 0.0, 0.0};
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dTheta, dPhi, dPsi;
  double cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  unsigned short iMarker, iVertex, nDim(3);
  unsigned long iPoint;

  double deltaT, q_n, qdot_n, qdot_nM1, q_nP1;
  double alpha_n, alphadot_n, alphadot_nM1, alpha_nP1;
  double alpha0, alpha1;
  double disp, dAlpha;
  double varCoordNorm2(0.0);

  //--- Get the current center of rotation (can vary at each iteration) ---
  Center[0] = structure->GetCenterOfRotation_x();
  Center[1] = structure->GetCenterOfRotation_y();
  Center[2] = structure->GetCenterOfRotation_z();

  //--- Get the center of rotation from previous time step ---
  Center_n[0] = structure->GetCenterOfRotation_n_x();
  Center_n[1] = structure->GetCenterOfRotation_n_y();
  Center_n[2] = structure->GetCenterOfRotation_n_z();

  if(!prediction){
    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -( (*(integrator->GetDisp()))[1] - (*(integrator->GetDisp_n()))[1]);
      newCenter[0] = Center[0];
      newCenter[1] = -(*(integrator->GetDisp()))[0];
      newCenter[2] = Center[2];
      disp = newCenter[1] - Center_n[1];
    }
    else if(config->GetStructType() == "SPRING_HOR"){
      dPsi = 0.0;
      newCenter[0] = (*(integrator->GetDisp()))[0];
      newCenter[1] = Center[1];
      newCenter[2] = Center[2];
    }
  }
  else{
    deltaT = config->GetDeltaT();
    alpha0 = 1.0;
    alpha1 = 0.5; //Second order prediction

    q_n = (*(integrator->GetDisp()))[0];
    qdot_n = (*(integrator->GetVel()))[0];
    qdot_nM1 = (*(integrator->GetVel_n()))[0];
    q_nP1 = q_n + alpha0*deltaT*qdot_n + alpha1*deltaT*(qdot_n - qdot_nM1);
    disp = q_nP1-q_n;

    if(structure->GetnDof() == 2){
      alpha_n = (*(integrator->GetDisp()))[1];
      alphadot_n = (*(integrator->GetVel()))[1];
      alphadot_nM1 = (*(integrator->GetVel_n()))[1];
      alpha_nP1 = alpha_n + alpha0*deltaT*alphadot_n + alpha1*deltaT*(alphadot_n - alphadot_nM1);
      dAlpha = alpha_nP1-alpha_n;
    }
    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -dAlpha;
      newCenter[0] = Center[0];
      newCenter[1] = -q_nP1;
      newCenter[2] = Center[2];
    }
    else{
      dPsi = 0.0;
    }
  }

  cosTheta = cos(dTheta);  cosPhi = cos(dPhi);  cosPsi = cos(dPsi);
  sinTheta = sin(dTheta);  sinPhi = sin(dPhi);  sinPsi = sin(dPsi);

  //--- Compute the rotation matrix. The implicit ordering is rotation about the x-axis, y-axis, then z-axis. ---

  rotMatrix[0][0] = cosPhi*cosPsi;
  rotMatrix[1][0] = cosPhi*sinPsi;
  rotMatrix[2][0] = -sinPhi;

  rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
  rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
  rotMatrix[2][1] = sinTheta*cosPhi;

  rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
  rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
  rotMatrix[2][2] = cosTheta*cosPhi;

  for(iMarker = 0; iMarker < geometry->GetnMarkers(); iMarker++){
    if (geometry->markersMoving[iMarker] == true){
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex];
        Coord = geometry->node[iPoint]->GetCoord();
        Coord_n = geometry->node[iPoint]->GetCoord_n();

        if(config->GetUnsteady() == "YES"){
          for (int iDim=0; iDim < nDim; iDim++){
            r[iDim] = Coord_n[iDim] - Center_n[iDim];
          }
        }
        else{
          for (int iDim=0; iDim < nDim; iDim++){
            r[iDim] = Coord[iDim] - Center[iDim];
          }
        }

        rotCoord[0] = rotMatrix[0][0]*r[0]
              + rotMatrix[0][1]*r[1]
              + rotMatrix[0][2]*r[2];

        rotCoord[1] = rotMatrix[1][0]*r[0]
              + rotMatrix[1][1]*r[1]
              + rotMatrix[1][2]*r[2];

        rotCoord[2] = rotMatrix[2][0]*r[0]
              + rotMatrix[2][1]*r[1]
              + rotMatrix[2][2]*r[2];

        for(int iDim=0; iDim < nDim; iDim++){
                newCoord[iDim] = newCenter[iDim] + rotCoord[iDim];
                varCoord[iDim] = newCoord[iDim] - Coord[iDim];
        }

        varCoordNorm2 += varCoord[0]*varCoord[0] + varCoord[1]*varCoord[1] + varCoord[2]*varCoord[2];

        //--- Apply change of coordinates to the node on the moving interface ---
        geometry->node[iPoint]->SetCoord(newCoord);

        //--- At initialisation, propagate the initial position of the inteface in the past ---
        if(initialize) geometry->node[iPoint]->SetCoord_n(newCoord);
      }
    }
  }

  varCoordNorm = sqrt(varCoordNorm2);

  //--- Update the position of the center of rotation ---
  structure->SetCenterOfRotation_X(newCenter[0]);
  structure->SetCenterOfRotation_Y(newCenter[1]);
  structure->SetCenterOfRotation_Z(newCenter[2]);
}
*/

void NativeSolidSolver::computeInterfacePosVel(bool initialize){

    double *Coord, *Coord_n, newCoord[3], newVel[3], Center[3], Center_n[3], newCenter[3], centerVel[3], rotCoord[3], r[3];
    double varCoord[3] = {0.0, 0.0, 0.0};
    double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    double dTheta, dPhi, dPsi;
    double psidot;
    double cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
    unsigned short iMarker, iVertex, nDim(3);
    unsigned long iPoint;
    double varCoordNorm2(0.0);

    /*--- Get the current center of rotation (can vary at each iteration) ---*/
    Center[0] = structure->GetCenterOfRotation_x();
    Center[1] = structure->GetCenterOfRotation_y();
    Center[2] = structure->GetCenterOfRotation_z();

    /*--- Get the center of rotation from previous time step ---*/
    Center_n[0] = structure->GetCenterOfRotation_n_x();
    Center_n[1] = structure->GetCenterOfRotation_n_y();
    Center_n[2] = structure->GetCenterOfRotation_n_z();

    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -( (integrator->GetSolver()->GetDisp())[1] - (integrator->GetSolver()->GetDisp_n())[1] );
      psidot = (integrator->GetSolver()->GetVel())[1];
      newCenter[0] = Center[0];
      newCenter[1] = -(integrator->GetSolver()->GetDisp())[0];
      newCenter[2] = Center[2];
      centerVel[0] = 0.0;
      centerVel[1] = -(integrator->GetSolver()->GetVel())[0];
      centerVel[2] = 0.0;
    }
    else if(config->GetStructType() == "SPRING_HOR"){
      dPsi = 0.0;
      psidot = 0.0;
      newCenter[0] = (integrator->GetSolver()->GetDisp())[0];
      newCenter[1] = Center[1];
      newCenter[2] = Center[2];
      centerVel[0] = (integrator->GetSolver()->GetVel())[0];
      centerVel[1] = 0.0;
      centerVel[2] = 0.0;
    }
    else if(config->GetStructType() == "SPRING_VER"){
      dPsi = 0.0;
      psidot = 0.0;
      newCenter[0] = Center[0];
      newCenter[1] = (integrator->GetSolver()->GetDisp())[0];
      newCenter[2] = Center[2];
      centerVel[0] = 0.0;
      centerVel[1] = (integrator->GetSolver()->GetVel())[0];
      centerVel[2] = 0.0;
    }
    else { }

    cosTheta = cos(dTheta);  cosPhi = cos(dPhi);  cosPsi = cos(dPsi);
    sinTheta = sin(dTheta);  sinPhi = sin(dPhi);  sinPsi = sin(dPsi);

    /*--- Compute the rotation matrix. The implicit
    ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

    rotMatrix[0][0] = cosPhi*cosPsi;
    rotMatrix[1][0] = cosPhi*sinPsi;
    rotMatrix[2][0] = -sinPhi;

    rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotMatrix[2][1] = sinTheta*cosPhi;

    rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotMatrix[2][2] = cosTheta*cosPhi;


    for(iMarker = 0; iMarker < geometry->GetnMarkers(); iMarker++){
      if (geometry->markersMoving[iMarker] == true){
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){

          iPoint = geometry->vertex[iMarker][iVertex];
          Coord = geometry->node[iPoint]->GetCoord();
          Coord_n = geometry->node[iPoint]->GetCoord_n();

          if(config->GetUnsteady() == "YES"){
            for (int iDim=0; iDim < nDim; iDim++){
              r[iDim] = Coord_n[iDim] - Center_n[iDim];
            }
          }
          else{
            for (int iDim=0; iDim < nDim; iDim++){
              r[iDim] = Coord[iDim] - Center[iDim];
            }
          }

          rotCoord[0] = rotMatrix[0][0]*r[0]
                + rotMatrix[0][1]*r[1]
                + rotMatrix[0][2]*r[2];

          rotCoord[1] = rotMatrix[1][0]*r[0]
                + rotMatrix[1][1]*r[1]
                + rotMatrix[1][2]*r[2];

          rotCoord[2] = rotMatrix[2][0]*r[0]
                + rotMatrix[2][1]*r[1]
                + rotMatrix[2][2]*r[2];

          for(int iDim=0; iDim < nDim; iDim++){
                  newCoord[iDim] = newCenter[iDim] + rotCoord[iDim];
                  varCoord[iDim] = newCoord[iDim] - Coord[iDim];
          }

          newVel[0] = centerVel[0] + psidot*(newCoord[1]-newCenter[1]);
          newVel[1] = centerVel[1] - psidot*(newCoord[0]-newCenter[0]);
          newVel[2] = centerVel[2] + 0.0;

          varCoordNorm2 += varCoord[0]*varCoord[0] + varCoord[1]*varCoord[1] + varCoord[2]*varCoord[2];

          /*--- Apply change of coordinates to the node on the moving interface ---*/
          geometry->node[iPoint]->SetCoord(newCoord);
          geometry->node[iPoint]->SetVel(newVel);

          /*--- At initialisation, propagate the initial position of the inteface in the past ---*/
          if(initialize && config->GetUnsteady() != "HARMONIC"){
              geometry->node[iPoint]->SetCoord_n(newCoord);
              geometry->node[iPoint]->SetVel_n(newVel);
          }
        }
      }
    }

    varCoordNorm = sqrt(varCoordNorm2);

    /*--- Update the position of the center of rotation ---*/
    structure->SetCenterOfRotation_X(newCenter[0]);
    structure->SetCenterOfRotation_Y(newCenter[1]);
    structure->SetCenterOfRotation_Z(newCenter[2]);

}

void NativeSolidSolver::computeInterfacePosVel(bool initialize, unsigned int instance){

    double *Coord, *Coord_n, newCoord[3], newVel[3], Center[3], Center_n[3], newCenter[3], centerVel[3], rotCoord[3], r[3];
    double varCoord[3] = {0.0, 0.0, 0.0};
    double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    double dTheta, dPhi, dPsi;
    double psidot;
    double cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
    unsigned short iMarker, iVertex, nDim(3);
    unsigned long iPoint;
    unsigned int nInstances = 2*config->GetNumberHarmonics()+1;
    double varCoordNorm2(0.0);

    if(config->GetUnsteady() != "HARMONIC"){
      instance = 0;
      nInstances = 1;
    }
    else if (instance>=nInstances){
      instance = nInstances-1;
    }
    

    /*--- Get the current center of rotation (can vary at each iteration) ---*/
    Center[0] = structure->GetCenterOfRotation_x();
    Center[1] = structure->GetCenterOfRotation_y();
    Center[2] = structure->GetCenterOfRotation_z();

    /*--- Get the center of rotation from previous time step ---*/
    Center_n[0] = structure->GetCenterOfRotation_n_x();
    Center_n[1] = structure->GetCenterOfRotation_n_y();
    Center_n[2] = structure->GetCenterOfRotation_n_z();

    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -( (integrator->GetSolver()->GetDisp())[instance+nInstances] - (integrator->GetSolver()->GetDisp_n())[instance+nInstances] );
      psidot = (integrator->GetSolver()->GetVel())[instance+nInstances];
      newCenter[0] = Center[0];
      newCenter[1] = -(integrator->GetSolver()->GetDisp())[instance];
      newCenter[2] = Center[2];
      centerVel[0] = 0.0;
      centerVel[1] = -(integrator->GetSolver()->GetVel())[instance];
      centerVel[2] = 0.0;
    }
    else if(config->GetStructType() == "SPRING_HOR"){
      dPsi = 0.0;
      psidot = 0.0;
      newCenter[0] = (integrator->GetSolver()->GetDisp())[instance];
      newCenter[1] = Center[1];
      newCenter[2] = Center[2];
      centerVel[0] = (integrator->GetSolver()->GetVel())[instance];
      centerVel[1] = 0.0;
      centerVel[2] = 0.0;
    }
    else if(config->GetStructType() == "SPRING_VER"){
      dPsi = 0.0;
      psidot = 0.0;
      newCenter[0] = Center[0];
      newCenter[1] = (integrator->GetSolver()->GetDisp())[instance];
      newCenter[2] = Center[2];
      centerVel[0] = 0.0;
      centerVel[1] = (integrator->GetSolver()->GetVel())[instance];
      centerVel[2] = 0.0;
    }
    else { }

    cosTheta = cos(dTheta);  cosPhi = cos(dPhi);  cosPsi = cos(dPsi);
    sinTheta = sin(dTheta);  sinPhi = sin(dPhi);  sinPsi = sin(dPsi);

    /*--- Compute the rotation matrix. The implicit
    ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

    rotMatrix[0][0] = cosPhi*cosPsi;
    rotMatrix[1][0] = cosPhi*sinPsi;
    rotMatrix[2][0] = -sinPhi;

    rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotMatrix[2][1] = sinTheta*cosPhi;

    rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotMatrix[2][2] = cosTheta*cosPhi;


    for(iMarker = 0; iMarker < geometry->GetnMarkers(); iMarker++){
      if (geometry->markersMoving[iMarker] == true){
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){

          iPoint = geometry->vertex[iMarker][iVertex];
          Coord = geometry->node[iPoint]->GetCoord();
          Coord_n = geometry->node[iPoint]->GetCoord_n();

          if(config->GetUnsteady() == "YES" || config->GetUnsteady() == "HARMONIC"){
            for (int iDim=0; iDim < nDim; iDim++){
              r[iDim] = Coord_n[iDim] - Center_n[iDim];
            }
          }
          else{
            for (int iDim=0; iDim < nDim; iDim++){
              r[iDim] = Coord[iDim] - Center[iDim];
            }
          }

          rotCoord[0] = rotMatrix[0][0]*r[0]
                + rotMatrix[0][1]*r[1]
                + rotMatrix[0][2]*r[2];

          rotCoord[1] = rotMatrix[1][0]*r[0]
                + rotMatrix[1][1]*r[1]
                + rotMatrix[1][2]*r[2];

          rotCoord[2] = rotMatrix[2][0]*r[0]
                + rotMatrix[2][1]*r[1]
                + rotMatrix[2][2]*r[2];

          for(int iDim=0; iDim < nDim; iDim++){
                  newCoord[iDim] = newCenter[iDim] + rotCoord[iDim];
                  varCoord[iDim] = newCoord[iDim] - Coord[iDim];
          }

          newVel[0] = centerVel[0] + psidot*(newCoord[1]-newCenter[1]);
          newVel[1] = centerVel[1] - psidot*(newCoord[0]-newCenter[0]);
          newVel[2] = centerVel[2] + 0.0;

          varCoordNorm2 += varCoord[0]*varCoord[0] + varCoord[1]*varCoord[1] + varCoord[2]*varCoord[2];

          /*--- Apply change of coordinates to the node on the moving interface ---*/
          geometry->node[iPoint]->SetCoord(newCoord);
          geometry->node[iPoint]->SetVel(newVel);

          /*--- At initialisation, propagate the initial position of the inteface in the past ---*/
          if(initialize){
              geometry->node[iPoint]->SetCoord_n(newCoord);
              geometry->node[iPoint]->SetVel_n(newVel);
          }
        }
      }
    }

    varCoordNorm = sqrt(varCoordNorm2);

    /*--- Update the position of the center of rotation ---*/
    structure->SetCenterOfRotation_X(newCenter[0]);
    structure->SetCenterOfRotation_Y(newCenter[1]);
    structure->SetCenterOfRotation_Z(newCenter[2]);

}

void NativeSolidSolver::setInitialDisplacements(){
  //mapRigidBodyMotion(false, true);
  computeInterfacePosVel(true);
}


void NativeSolidSolver::staticComputation(){

  integrator->StaticIteration(structure);

  //mapRigidBodyMotion(false,false);
  computeInterfacePosVel(false);
}

void NativeSolidSolver::harmonicComputation(){

  integrator->HarmonicIteration(config, structure);

  //mapRigidBodyMotion(false,false);
  computeInterfacePosVel(false);
}


void NativeSolidSolver::writeSolution(double currentTime, double lastTime, double currentFSIIter, unsigned long ExtIter, unsigned long NbExtIter){

  int rank = MASTER_NODE;
  double DeltaT;
  string restartFileName;
  unsigned long DeltaIter = config->GetDeltaIterWrite();

  DeltaT = currentTime-lastTime;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if(rank == MASTER_NODE){
  if(structure->GetnDof() == 1){
    if(config->GetUnsteady() == "YES"){
      //if(currentTime == config->GetStartTime()) historyFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << "\"Acceleration variable\"" << endl;
      historyFile << currentTime << "\t" << (integrator->GetSolver()->GetDisp())[0] << "\t" << (integrator->GetSolver()->GetVel())[0] << "\t" << (integrator->GetSolver()->GetAcc())[0] << "\t" << (integrator->GetSolver()->GetAccVar())[0] << endl;
    }
    else{
      //if(currentFSIIter == config->GetStartTime()) historyFile << "\"FSI Iteration\"" << "\t" << "\"Displacement\"" << endl;
      historyFile << currentFSIIter << "\t" << (integrator->GetSolver()->GetDisp())[0] << endl;
    }
  }
  else if(structure->GetnDof() == 2){
    if(config->GetUnsteady() == "YES"){
      //if(currentTime == config->GetStartTime()) historyFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << "\t" << "\"Acceleration variable 1\"" << "\t" << "\"Acceleration variable 2\"" << endl;
      historyFile << currentTime << "\t" << (integrator->GetSolver()->GetDisp())[0] << "\t" << (integrator->GetSolver()->GetDisp())[1] << "\t" << (integrator->GetSolver()->GetVel())[0] << "\t" << (integrator->GetSolver()->GetVel())[1] << "\t" << (integrator->GetSolver()->GetAcc())[0] << "\t" << (integrator->GetSolver()->GetAcc())[1] << "\t" << (integrator->GetSolver()->GetAccVar())[0] << "\t" << (integrator->GetSolver()->GetAccVar())[1] << endl;
    }
    else{
      //if(currentFSIIter == config->GetStartTime()) historyFile << "\"FSI Iteration\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 1\"" << endl;
      historyFile << currentFSIIter << "\t" << (integrator->GetSolver()->GetDisp())[0] << "\t" << (integrator->GetSolver()->GetDisp())[1] << endl;
    }
  }
  }


  if(config->GetUnsteady() == "YES" || config->GetKindProblem() == "ADJOINT"){
  if ( ExtIter%DeltaIter == 0 || ExtIter == NbExtIter ){
    restartFileName = config->GetRestartFile();
    restartFile.open(restartFileName.c_str(), ios::out);
    if (structure->GetnDof() == 1){
      restartFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << "\"Acceleration variable\"" << endl;
      restartFile << currentTime-DeltaT << "\t" << (integrator->GetSolver()->GetDisp_n())[0] << "\t" << (integrator->GetSolver()->GetVel_n())[0] << "\t" << (integrator->GetSolver()->GetAcc_n())[0] << "\t" << (integrator->GetSolver()->GetAccVar_n())[0] << endl;
      restartFile << currentTime << "\t" << (integrator->GetSolver()->GetDisp())[0] << "\t" << (integrator->GetSolver()->GetVel())[0] << "\t" << (integrator->GetSolver()->GetAcc())[0] << "\t" << (integrator->GetSolver()->GetAccVar())[0] << endl;
    }
    else if (structure->GetnDof() == 2){
      restartFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << "\t" << "\"Acceleration variable 1\"" << "\t" << "\"Acceleration variable 2\"" << endl;
      restartFile << currentTime-DeltaT << "\t" << (integrator->GetSolver()->GetDisp_n())[0] << "\t" << (integrator->GetSolver()->GetDisp_n())[1] << "\t" << (integrator->GetSolver()->GetVel_n())[0] << "\t" << (integrator->GetSolver()->GetVel_n())[1] << "\t" << (integrator->GetSolver()->GetAcc_n())[0] << "\t" << (integrator->GetSolver()->GetAcc_n())[1] << "\t" << (integrator->GetSolver()->GetAccVar_n())[0] << "\t" << (integrator->GetSolver()->GetAccVar_n())[1] << endl;
      restartFile << currentTime << "\t" << (integrator->GetSolver()->GetDisp())[0] << "\t" << (integrator->GetSolver()->GetDisp())[1] << "\t" << (integrator->GetSolver()->GetVel())[0] << "\t" << (integrator->GetSolver()->GetVel())[1] << "\t" << (integrator->GetSolver()->GetAcc())[0] << "\t" << (integrator->GetSolver()->GetAcc())[1] << "\t" << (integrator->GetSolver()->GetAccVar())[0] << "\t" << (integrator->GetSolver()->GetAccVar())[1] << endl;
    }
    restartFile.close();
  }
  }

}

void NativeSolidSolver::saveSolution(){

    int rank = MASTER_NODE;
    double DeltaT;
    //string restartFileName;
    //unsigned long DeltaIter = config->GetDeltaIterWrite();
    unsigned long ExtIter = integrator->GetExtIter();

    DeltaT = config->GetDeltaT();

  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

    if(rank == MASTER_NODE){
      if(structure->GetnDof() == 1){
        if(config->GetUnsteady() == "YES"){
          historyFile << fixed
                      << setw(10) << DeltaT
                      << setw(10) << ExtIter
                      << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                      << setw(15) << (integrator->GetSolver()->GetVel())[0]
                      << setw(15) << (integrator->GetSolver()->GetAcc())[0] << endl;
        }
        else{
          historyFile << fixed
                      << setw(10) << DeltaT
                      << setw(10) << ExtIter
                      << setw(15) << (integrator->GetSolver()->GetDisp())[0] << endl;
        }
      }
      else if(structure->GetnDof() == 2){
        if(config->GetUnsteady() == "YES"){
          historyFile << fixed
                      << setw(10) << DeltaT
                      << setw(10) << ExtIter
                      << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                      << setw(15) << (integrator->GetSolver()->GetDisp())[1]
                      << setw(15) << (integrator->GetSolver()->GetVel())[0]
                      << setw(15) << (integrator->GetSolver()->GetVel())[1]
                      << setw(15) << (integrator->GetSolver()->GetAcc())[0]
                      << setw(15) << (integrator->GetSolver()->GetAcc())[1] << endl;
        }
        else{
          historyFile << fixed
                      << setw(10) << DeltaT
                      << setw(10) << ExtIter
                      << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                      << setw(15) << (integrator->GetSolver()->GetDisp())[1] << endl;
        }
      }
    }
  output->WriteRestart(integrator, structure, config);
}

void NativeSolidSolver::writeSolution(double time, int FSIter){

    int rank = MASTER_NODE;

  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

    if(rank == MASTER_NODE){
    if(structure->GetnDof() == 1){
      if(config->GetUnsteady() == "YES"){
        historyFile2 << fixed
                     << setw(10) << time
                     << setw(10) << FSIter
                     << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                     << setw(15) << (integrator->GetSolver()->GetVel())[0]
                     << setw(15) << (integrator->GetSolver()->GetAcc())[0] << endl;
      }
      else if(config->GetUnsteady() == "HARMONIC"){
        double currentTime = 0.;
        unsigned short nInst = 2*getNumberHarmonics()+1;
        for (unsigned short iInst = 0; iInst < nInst ; iInst++)
        {
          historyFile2 << fixed
                       << setw(10) << currentTime
                       << setw(10) << FSIter
                       << setw(15) << (integrator->GetSolver()->GetDisp())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetVel())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetAcc())[iInst] << endl;
          currentTime += config->GetStopTime()/nInst;
        }
      }
      else{
        historyFile2 << fixed
                     << setw(10) << FSIter
                     << setw(15) << (integrator->GetSolver()->GetDisp())[0] << endl;
      }
    }
    else if(structure->GetnDof() == 2){
      if(config->GetUnsteady() == "YES"){
        historyFile2 << fixed
                     << setw(10) << time
                     << setw(10) << FSIter
                     << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                     << setw(15) << (integrator->GetSolver()->GetDisp())[1]
                     << setw(15) << (integrator->GetSolver()->GetVel())[0]
                     << setw(15) << (integrator->GetSolver()->GetVel())[1]
                     << setw(15) << (integrator->GetSolver()->GetAcc())[0]
                     << setw(15) << (integrator->GetSolver()->GetAcc())[1] << endl;
      }
      else if(config->GetUnsteady() == "HARMONIC"){
        double currentTime = 0.;
        double Energy(0);
        unsigned short nInst = 2*getNumberHarmonics()+1;
        for (unsigned short iInst = 0; iInst < nInst ; iInst++)
        {
          double h, alpha;
          h = (integrator->GetSolver()->GetDisp())[iInst];
          alpha = (integrator->GetSolver()->GetDisp())[iInst+nInst];
          Energy += .5*(structure->Get_Kh())*h*h;
          Energy += .5*(structure->Get_Ka())*alpha*alpha;
          historyFile2 << fixed
                       << setw(10) << currentTime
                       << setw(10) << FSIter
                       << setw(15) << (integrator->GetSolver()->GetDisp())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetDisp())[iInst+nInst]
                       << setw(15) << (integrator->GetSolver()->GetVel())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetVel())[iInst+nInst]
                       << setw(15) << (integrator->GetSolver()->GetAcc())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetAcc())[iInst+nInst]
                       << setw(15) << Energy
                       << setw(15) << (integrator->GetSolver()->GetOmega()) 
                       << setw(15) << (integrator->GetSolver()->GetDeltaOmega()) << endl;
          currentTime += config->GetStopTime()/nInst;
        }
      }
      else{
        historyFile2 << fixed
                     << setw(10) << FSIter
                     << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                     << setw(15) << (integrator->GetSolver()->GetDisp())[1] << endl;
      }
    }
    }
}

void NativeSolidSolver::writeAdjointSolution(double time, int FSIter){

    int rank = MASTER_NODE;
    double dampingDer = 0.0;

  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

    if(rank == MASTER_NODE){
    if(structure->GetnDof() == 1){
      if(config->GetUnsteady() == "YES"){
        historyFile2 << fixed
                     << setw(10) << time
                     << setw(10) << FSIter
                     << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                     << setw(15) << (integrator->GetSolver()->GetVel())[0]
                     << setw(15) << (integrator->GetSolver()->GetAcc())[0] << endl;
      }
      else if(config->GetUnsteady() == "HARMONIC"){
        double currentTime = 0.;
        double Energy(0);
        double h(0), v(0), a(0), dEdkh(0), dRdkh(0), dEdm(0), dRdm(0);
        unsigned short nInst = 2*getNumberHarmonics()+1;

        for (unsigned short iInst = 0; iInst < nInst ; iInst++){
          h = (integrator->GetSolver()->GetDisp())[iInst];
          v = (integrator->GetSolver()->GetVel())[iInst];
          a = (integrator->GetSolver()->GetAcc())[iInst];
          dRdkh = h;
          Energy += .5*(structure->Get_Kh())*h*h;
          dEdkh += .5*h*h-(integrator->GetSolver()->GetAdjLoads())[iInst]*dRdkh;
          dRdm = a;
          dEdm += .5*v*v-(integrator->GetSolver()->GetAdjLoads())[iInst]*dRdm;
        }
        dampingDer = getPlungeDampingDerivative();
        for (unsigned short iInst = 0; iInst < nInst ; iInst++)
        {
          historyFile2 << fixed
                       << setw(10) << currentTime
                       << setw(10) << FSIter
                       << setw(15) << (integrator->GetSolver()->GetDisp())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetVel())[iInst]
                       << setw(15) << (integrator->GetSolver()->GetAcc())[iInst] << endl;
          currentTime += config->GetStopTime()/nInst;
        }
      }
      else{
        double Energy(0);
        double h(0), dEdkh(0), dRdkh(0);
        h = (integrator->GetSolver()->GetDisp())[0];
        dRdkh = h;
        Energy += .5*(structure->Get_Kh())*h*h;
        dEdkh = .5*h*h-(integrator->GetSolver()->GetAdjLoads())[0]*dRdkh;
        historyFile2 << fixed
                     << setw(10) << FSIter
                     << setw(15) << Energy
                     << setw(15) << dEdkh << endl;
      }
    }
    else if(structure->GetnDof() == 2){
      if(config->GetUnsteady() == "YES"){
        historyFile2 << fixed
                     << setw(10) << time
                     << setw(10) << FSIter
                     << setw(15) << (integrator->GetSolver()->GetDisp())[0]
                     << setw(15) << (integrator->GetSolver()->GetDisp())[1]
                     << setw(15) << (integrator->GetSolver()->GetVel())[0]
                     << setw(15) << (integrator->GetSolver()->GetVel())[1]
                     << setw(15) << (integrator->GetSolver()->GetAcc())[0]
                     << setw(15) << (integrator->GetSolver()->GetAcc())[1] << endl;
      }
      else if(config->GetUnsteady() == "HARMONIC"){
        double Energy(0), dJdkh(0), dJdka(0);
        double h, alpha;
        double currentTime = 0.;
        unsigned short nInst = 2*getNumberHarmonics()+1;
        for (unsigned short iInst = 0; iInst < nInst ; iInst++)
        {
          h = (integrator->GetSolver()->GetDisp())[iInst];
          alpha = (integrator->GetSolver()->GetDisp())[iInst+nInst];
          Energy += .5*(structure->Get_Kh())*h*h;
          Energy += .5*(structure->Get_Ka())*alpha*alpha;
          if (config->GetObjFunction() == "ENERGY")
          {
            dJdkh += .5*h*h;
            dJdka += .5*alpha*alpha;
          }
          dampingDer = getPlungeDampingDerivative();
          historyFile2 << fixed
                       << setw(10) << currentTime
                       << setw(10) << FSIter
                       << setw(15) << Energy
                       << setw(15) << integrator->GetSolver()->GetStiffnessDerivative(0)+dJdkh
                       << setw(15) << integrator->GetSolver()->GetStiffnessDerivative(1)+dJdka
                       << setw(15) << dampingDer
                       << setw(15) << integrator->GetSolver()->GetDampingDerivative(1)
                       << setw(15) << integrator->GetSolver()->GetMassDerivative(0)+integrator->GetSolver()->GetImbalanceDerivative()*structure->Get_S()/structure->Get_m()
                       << setw(15) << integrator->GetSolver()->GetMassDerivative(1)
                       << setw(15) << integrator->GetSolver()->GetImbalanceDerivative()/structure->Get_m() << endl;
          historyFile << scientific
                       << setw(10) << currentTime
                       << setw(10) << FSIter
                       << setw(20) << Energy
                       << setw(20) << integrator->GetSolver()->GetAdjLoads()[iInst]
                       << setw(20) << integrator->GetSolver()->GetAdjLoads()[iInst+nInst]<< endl;
          currentTime += config->GetStopTime()/nInst;
        }
      }
      else{
        double Energy(0);
        double h(0), dEdkh(0), dRdkh(0);
        double alpha(0), dEdka(0), dRdka(0);
        h = (integrator->GetSolver()->GetDisp())[0];
        dRdkh = h;
        dEdkh = .5*h*h-(integrator->GetSolver()->GetAdjLoads())[0]*dRdkh;
        alpha = (integrator->GetSolver()->GetDisp())[1];
        dRdka = alpha;
        dEdka = .5*alpha*alpha-(integrator->GetSolver()->GetAdjLoads())[1]*dRdka;
        Energy += .5*(structure->Get_Kh())*h*h;
        Energy += .5*(structure->Get_Ka())*alpha*alpha;
        historyFile2 << fixed
                     << setw(10) << FSIter
                     << setw(15) << Energy
                     << setw(15) << dEdkh
                     << setw(15) << dEdka << endl;
      }
    }
  }
  output->WriteAdjointOutput(dampingDer);
}

void NativeSolidSolver::updateSolution(){

  if(config->GetUnsteady() == "YES"){
    integrator->UpdateSolution();
    geometry->UpdateGeometry();
    structure->SetCenterOfRotation_n_X(structure->GetCenterOfRotation_x());
    structure->SetCenterOfRotation_n_Y(structure->GetCenterOfRotation_y());
    structure->SetCenterOfRotation_n_Z(structure->GetCenterOfRotation_z());
  }
  else
    q_uM1 = (integrator->GetSolver()->GetDisp());
}

/*
 * void NativeSolidSolver::updateGeometry(){

  if(config->GetUnsteady() == "YES"){
    geometry->UpdateGeometry();
    structure->SetCenterOfRotation_n_X(structure->GetCenterOfRotation_x());
    structure->SetCenterOfRotation_n_Y(structure->GetCenterOfRotation_y());
    structure->SetCenterOfRotation_n_Z(structure->GetCenterOfRotation_z());
  }
}*/

/*
void NativeSolidSolver::displacementPredictor(){

  mapRigidBodyMotion(true, false);

}*/

unsigned short NativeSolidSolver::getFSIMarkerID(){

  unsigned short iMarker(0), IDtoSend;


  while(iMarker < geometry->GetnMarkers()){
      if(geometry->GetMarkersMoving(iMarker)){
          IDtoSend =  iMarker;
          break;
      }
      iMarker++;
  }
  return IDtoSend;
}

unsigned long NativeSolidSolver::getNumberOfSolidInterfaceNodes(unsigned short iMarker){
  return nSolidInterfaceVertex;
}

unsigned int NativeSolidSolver::getInterfaceNodeGlobalIndex(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;

  iPoint = geometry->vertex[iMarker][iVertex];

  return iPoint;
}

double NativeSolidSolver::getInterfaceNodePosX(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();

    return Coord[0];
}

double NativeSolidSolver::getInterfaceNodePosY(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();

    return Coord[1];
}

double NativeSolidSolver::getInterfaceNodePosZ(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();

    return 0.0; //3D is not really implemented in this solver...
}

double NativeSolidSolver::getInterfaceNodePosX0(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord0();

    return Coord[0];
}

double NativeSolidSolver::getInterfaceNodePosY0(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord0();

    return Coord[1];
}

double NativeSolidSolver::getInterfaceNodePosZ0(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord0();

    return 0.0;
}

double NativeSolidSolver::getInterfaceNodeDispX(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord, *Coord0;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();
    Coord0 = geometry->node[iPoint]->GetCoord0();

    return Coord[0] - Coord0[0];
}

double NativeSolidSolver::getInterfaceNodeDispY(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord, *Coord0;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();
    Coord0 = geometry->node[iPoint]->GetCoord0();
    return Coord[1] - Coord0[1];
}

double NativeSolidSolver::getInterfaceNodeDispZ(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord, *Coord0;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();
    Coord0 = geometry->node[iPoint]->GetCoord0();

    return Coord[2] - Coord0[2];
}

double NativeSolidSolver::getInterfaceNodeVelX(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel = geometry->node[iPoint]->GetVel();

    return Vel[0];
}

double NativeSolidSolver::getInterfaceNodeVelY(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel = geometry->node[iPoint]->GetVel();

    return Vel[1];
}

double NativeSolidSolver::getInterfaceNodeVelZ(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel = geometry->node[iPoint]->GetVel();

    return Vel[2];
}

double NativeSolidSolver::getInterfaceNodeVelXNm1(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel_n;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel_n = geometry->node[iPoint]->GetVel_n();

    return Vel_n[0];
}

double NativeSolidSolver::getInterfaceNodeVelYNm1(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel_n;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel_n = geometry->node[iPoint]->GetVel_n();

    return Vel_n[1];
}

double NativeSolidSolver::getInterfaceNodeVelZNm1(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel_n;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel_n = geometry->node[iPoint]->GetVel_n();

    return Vel_n[2];
}

double NativeSolidSolver::getRotationCenterPosX(){

    return structure->GetCenterOfRotation_x();

}

double NativeSolidSolver::getRotationCenterPosY(){

    return structure->GetCenterOfRotation_y();
}

double NativeSolidSolver::getRotationCenterPosZ(){

    return structure->GetCenterOfRotation_z();
}

void NativeSolidSolver::setGeneralisedForce(unsigned int instance){

  unsigned short iVertex, iMarker;
  unsigned long iPoint, nPoint;
  double ForceX(0.0), ForceY(0.0), ForceZ(0.0);
  double* force;

  iMarker = getFSIMarkerID();
  nPoint = geometry->GetnPoint();

  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      iPoint = geometry->vertex[iMarker][iVertex];
      force = geometry->node[iPoint+nPoint*instance]->GetForce();
      ForceX += force[0];
      ForceY += force[1];
      ForceZ += force[2];
  }

  if(config->GetStructType() == "SPRING_HOR"){
    (integrator->GetSolver()->GetLoads())[instance] = ForceX;
  }
  else if(config->GetStructType() == "SPRING_VER"){
    (integrator->GetSolver()->GetLoads())[instance] = ForceY;
  }
  else if(config->GetStructType() == "AIRFOIL"){
    (integrator->GetSolver()->GetLoads())[instance] = -ForceY;
  }
  else{
    cerr << "Wrong structural type for applying global fluid loads !" << endl;
    throw(-1);
  }

}

void NativeSolidSolver::setGeneralisedForce(double Fx, double Fy){

  if(config->GetStructType() == "SPRING_HOR"){
    (integrator->GetSolver()->GetLoads())[0] = Fx;
  }
  else if(config->GetStructType() == "SPRING_VER"){
    (integrator->GetSolver()->GetLoads())[0] = Fy;
  }
  else if(config->GetStructType() == "AIRFOIL"){
    (integrator->GetSolver()->GetLoads())[0] = -Fy;
  }
  else{
    cerr << "Wrong structural type for applying global fluid loads !" << endl;
    throw(-1);
  }
}

void NativeSolidSolver::setGeneralisedMoment(unsigned int instance){

  unsigned short iVertex, iMarker;
  unsigned long iPoint, nPoint;
  double Moment(0.0), CenterX, CenterY, CenterZ;
  double dMoment(0.0); // Derivative of moment with respect to pitch
  double* Force;
  double* Coord;

  iMarker = getFSIMarkerID();
  nPoint = geometry->GetnPoint();

  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      iPoint = geometry->vertex[iMarker][iVertex];
      Force = geometry->node[iPoint+nPoint*instance]->GetForce();
      Coord = geometry->node[iPoint]->GetCoord();
      CenterX = getRotationCenterPosX();
      CenterY = getRotationCenterPosY();
      Moment += (Force[1]*(Coord[0]-CenterX) - Force[0]*(Coord[1]-CenterY));
      dMoment += Force[1]*(Coord[1]-CenterY);
      dMoment += Force[0]*(Coord[0]-CenterX); // TODO: Check signs
  }

  if(config->GetStructType() == "AIRFOIL"){
    (integrator->GetSolver()->GetLoads())[2*config->GetNumberHarmonics()+1+instance] = -Moment;
    integrator->GetSolver()->SetLoadGradient(instance, dMoment); // This is added to the Jacobian, so it's positive
  }
  else if(config->GetStructType() == "SPRING_VER"){}
  else if(config->GetStructType() == "SPRING_HOR"){}
  else{
    cerr << "Wrong structural type for applying global fluid loads !" << endl;
    throw(-1);
  }

}

void NativeSolidSolver::setGeneralisedMoment(double M){

  if(config->GetStructType() == "AIRFOIL"){
    (integrator->GetSolver()->GetLoads())[1] = -M;
  }
  else if(config->GetStructType() == "SPRING_VER"){}
  else if(config->GetStructType() == "SPRING_HOR"){}
  else{
    cerr << "Wrong structural type for applying global fluid loads !" << endl;
    throw(-1);
  }
}

void NativeSolidSolver::setTotalAdjointDisplacement(unsigned int instance){

  unsigned short iVertex, iMarker;
  unsigned long iPoint;
  double AdjDispX(0.0), AdjDispY(0.0), AdjDispZ(0.0);
  double AdjVelX(0.0), AdjVelY(0.0), AdjVelZ(0.0);
  double CenterX, CenterY, CenterZ;
  double* ExtAdjointDisp;
  double* Coord;
  double dPsi(0.0), sinPsi(0.0), cosPsi(1.0);
  unsigned int offset = instance+2*config->GetNumberHarmonics()+1;
  unsigned int offset_vel = offset+2*config->GetNumberHarmonics()+1;
  unsigned long nPoint = geometry->GetnPoint();
  
  if (config->GetStructType() == "AIRFOIL")
  {
    dPsi = (integrator->GetSolver()->GetDisp())[offset]; // dPsi = -alpha in the direct part but it gives me a headache
    cout << "dPsi: " << dPsi << endl;
    sinPsi = sin(dPsi);
    cosPsi = cos(dPsi);
  }
  

  iMarker = getFSIMarkerID();
  CenterX = getRotationCenterPosX();
  CenterY = getRotationCenterPosY();
  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      iPoint = geometry->vertex[iMarker][iVertex];
      ExtAdjointDisp = geometry->node[iPoint+nPoint*instance]->GetDisplacementAdjoint();
      AdjDispX += ExtAdjointDisp[0]; // dx/dx*dJ/dx
      AdjDispY += ExtAdjointDisp[1]; // dy/dy*dJ/dy
      Coord = geometry->node[iPoint]->GetCoord();
      AdjDispZ += ExtAdjointDisp[0]*(Coord[1]-CenterY); // rotation <- z, this is dx/dPsi*dJ/dx AROUND THE DEFORMED SHAPE. Therefore dPsi is infinitesimally small!
      AdjDispZ += ExtAdjointDisp[1]*(-1.0*(Coord[0]-CenterX)); // rotation <- z, this is dy/dPsi*dJ/dy
  }

  if (config->GetObjFunction() == "ENERGY")
  {
    AdjDispX += structure->Get_Kh()*integrator->GetSolver()->GetDisp()[instance];
    if (structure->GetnDof() == 2)
    {
      AdjDispY -= structure->Get_Kh()*integrator->GetSolver()->GetDisp()[instance]; // Negative because the sign of AdjDispY is negative for an airfoil
      AdjDispZ += structure->Get_Ka()*dPsi;
    }
    else
    {
      AdjDispY += structure->Get_Kh()*integrator->GetSolver()->GetDisp()[instance];
    }
    
  }
  else if (config->GetObjFunction() == "DISSIPATED_PLUNGE_POWER")
  {
    AdjVelX += 2*structure->Get_Ch()*integrator->GetSolver()->GetVel()[instance];
    AdjVelY += 2*structure->Get_Ch()*integrator->GetSolver()->GetVel()[instance];
  }
  
  


  cout << "Instance: " << instance << endl;
  cout << "AdjDispX: " << AdjDispX << endl;
  cout << "AdjDispY: " << AdjDispY << endl;
  cout << "AdjDispZ: " << AdjDispZ << endl;

  if(config->GetStructType() == "AIRFOIL"){
    (integrator->GetSolver()->GetAdjDisps())[instance] = -AdjDispY;
    (integrator->GetSolver()->GetAdjDisps())[offset] = AdjDispZ;
    (integrator->GetSolver()->GetAdjDisps())[offset_vel] = AdjVelY;
  }
  else if(config->GetStructType() == "SPRING_VER"){
    (integrator->GetSolver()->GetAdjDisps())[instance] = AdjDispY;
    (integrator->GetSolver()->GetAdjDisps())[offset] = AdjVelY;
  }
  else if(config->GetStructType() == "SPRING_HOR"){
    (integrator->GetSolver()->GetAdjDisps())[instance] = AdjDispX;
    (integrator->GetSolver()->GetAdjDisps())[offset] = AdjVelX;
  }
  else{
    cerr << "Wrong structural type for applying global fluid loads !" << endl;
    throw(-1);
  }

}

void NativeSolidSolver::setFrequencyDerivative(double dJdw){
  integrator->GetSolver()->SetFrequencyDerivative(dJdw);
}

void NativeSolidSolver::setDesignVariableCentre(double x, unsigned long iDV)
{
  if (iDV < config->GetNumberDesignVariables())
  {
    posDV[iDV] = x;
  }
}

void NativeSolidSolver::setDesignVariableMagnitude(double mag, unsigned long iDV)
{
  if (iDV < config->GetNumberDesignVariables())
  {
    magDV[iDV] = mag;
  }
}

void NativeSolidSolver::setDesignVariableSide(double side, unsigned long iDV)
{
  if (iDV < config->GetNumberDesignVariables())
  {
    sideDV[iDV] = side;
  }
}

void NativeSolidSolver::applyDesignVariables()
{
  unsigned short iMarker, iVertex;
  unsigned long iPoint;
  double newCoord[3];
  double* Coord0;

  iMarker = getFSIMarkerID();
  if (config->GetDesignVariableKind() == "HICKS_HENNE"){
    for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      iPoint = geometry->vertex[iMarker][iVertex];
      Coord0 = geometry->node[iPoint]->GetCoord0();
      
      if (Coord0[0] > 1.0e-5) // To avoid problems in LE
      {
        double upper = copysign(1.0, Coord0[1]); // Differentiate between suction and pressure side. Will fail for some aerofoils...
        double xk = Coord0[0]/config->GetCord();
        newCoord[0] = Coord0[0];
        newCoord[1] = Coord0[1];
        newCoord[2] = Coord0[2];
        for (unsigned long iDV = 0; iDV < config->GetNumberDesignVariables(); iDV++){
          double side;
          switch ((int) sideDV[iDV]){
          case 1: // Bump in extrados
            side = (upper+1.)/2.;
            break;
          case -1: // Bump in intrados
            side = (upper-1.)/2.;
            break;
          default: // Symmetric bump
            side = upper;
            break;
          }
          
          double ek = log10(0.5)/log10(posDV[iDV]);
          newCoord[0] = newCoord[0];
          newCoord[1] = newCoord[1]+side*magDV[iDV]*pow(sin( M_PI * pow(xk, ek)), 3.0);
          newCoord[2] = newCoord[2];
        }
        geometry->node[iPoint]->SetCoord_n(newCoord); // Hacky. I am using Coord_n as the new base coordinate, with the DVs applied
      }
    }
  }
  
}

double NativeSolidSolver::getDesignVariableDerivative(unsigned long iDV)
{
  unsigned short iMarker, iVertex;
  unsigned long iPoint, nPoint;
  unsigned long nInst = 2*config->GetNumberHarmonics()+1;
  double dJdDV = 0.;
  double dX, dY, dZ, dPsi;
  double* Coord0;
  double* ExtAdjointDisp;
  double* Force;
  CVector AdjLoads = integrator->GetSolver()->GetAdjLoads();

  iMarker = getFSIMarkerID();
  nPoint = geometry->GetnPoint();
  if (iDV < config->GetNumberDesignVariables() && config->GetDesignVariableKind() == "HICKS_HENNE")
  {
    double ek = log10(0.5)/log10(posDV[iDV]);
    for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      dX = 0.;
      dY = 0.;
      dZ = 0.;
      iPoint = geometry->vertex[iMarker][iVertex];
      Coord0 = geometry->node[iPoint]->GetCoord0(); // We are deforming the original shape
      if (Coord0[0] > 1.0e-5) // To avoid problems in LE
      {
        double upper = copysign(1.0, Coord0[1]); // Differentiate between suction and pressure side. Will fail for some aerofoils...
        double xk = Coord0[0]/config->GetCord();
        double side;
        switch ((int) sideDV[iDV]){
        case 1: // Bump in extrados
          side = (upper+1.)/2.;
          break;
        case -1: // Bump in intrados
          side = (upper-1.)/2.;
          break;
        default: // Symmetric bump
          side = upper;
          break;
        }
        dX = 0.;
        dY = side*pow(sin( M_PI * pow(xk, ek)), 3.0);
        for (unsigned long iInst = 0; iInst < nInst; iInst++)
        {
          dPsi = (integrator->GetSolver()->GetDisp())[iInst+nInst];
          ExtAdjointDisp = geometry->node[iPoint+nPoint*iInst]->GetDisplacementAdjoint();
          Force = geometry->node[iPoint+nPoint*iInst]->GetForce();
          dJdDV += ExtAdjointDisp[0]*sin(dPsi)*dY; // SU2 dJ/dx*dx/dBump
          dJdDV += ExtAdjointDisp[1]*cos(dPsi)*dY; // SU2 dJ/dy*dy/dBump
          dJdDV -= AdjLoads[nInst+iInst]*Force[1]*sin(dPsi)*dY;
          dJdDV += AdjLoads[nInst+iInst]*Force[0]*cos(dPsi)*dY;
          /*Moment += (Force[1]*(Coord[0]-CenterX) - Force[0]*(Coord[1]-CenterY));
          dMoment += Force[1]*(Coord[1]-CenterY);
          dMoment += Force[0]*(Coord[0]-CenterX); // TODO: Check signs*/
        }
        
      }
    }

  }
  return dJdDV;
}

double NativeSolidSolver::getObjectiveFunction()
{
  double J = 0.0;
  unsigned short nInst = 2*getNumberHarmonics()+1;
  if (config->GetObjFunction() == "PITCH_AMPLITUDE")
  {
    J += integrator->GetSolver()->GetAmplitude();
  }
  else if (config->GetObjFunction() == "ENERGY")
  {
    for (unsigned short iInst = 0; iInst < nInst ; iInst++)
    {
      double h, alpha;
      h = (integrator->GetSolver()->GetDisp())[iInst];
      alpha = (integrator->GetSolver()->GetDisp())[iInst+nInst];
      J += .5*(structure->Get_Kh())*h*h;
      J += .5*(structure->Get_Ka())*alpha*alpha;
    }
  }
  else if (config->GetObjFunction() == "DISSIPATED_PLUNGE_POWER")
  {
    double hdot;
    for (unsigned short iInst = 0; iInst < nInst ; iInst++)
    {
      hdot = (integrator->GetSolver()->GetVel())[iInst];
      J += (structure->Get_Ch())*hdot*hdot;
    }
  }
  
  
  return J;
}

void NativeSolidSolver::applyload(unsigned short iVertex, unsigned int iInst, double Fx, double Fy, double Fz){

    unsigned short iMarker;
    unsigned long iPoint, nPoint;
    double Force[3];

    nPoint = geometry->GetnPoint();

    Force[0] = Fx;
    Force[1] = Fy;
    Force[2] = Fz;

    iMarker = getFSIMarkerID();
    iPoint = geometry->vertex[iMarker][iVertex];
    geometry->node[iPoint+nPoint*iInst]->SetForce(Force);

}

unsigned int NativeSolidSolver::getNumberHarmonics(){
    return config->GetNumberHarmonics();
}

void NativeSolidSolver::applyDisplacementAdjoint(unsigned short iVertex, unsigned int iInst, double Dx, double Dy, double Dz){
  if (config->GetKindProblem() == "ADJOINT")
  {
    unsigned short iMarker;
    unsigned long iPoint, nPoint;
    double DispAdj[3];

    DispAdj[0] = Dx;
    DispAdj[1] = Dy;
    DispAdj[2] = Dz;

    iMarker = getFSIMarkerID();
    iPoint = geometry->vertex[iMarker][iVertex];
    nPoint = geometry->GetnPoint();
    geometry->node[iPoint+iInst*nPoint]->SetDisplacementAdjoint(DispAdj);
  }
  
}

double NativeSolidSolver::getLoadDerivativeX(unsigned short iVertex){
  if (config->GetKindProblem() == "ADJOINT")
  {
    unsigned short iMarker;
    unsigned long iPoint;
    double *LoadDerivative;
    iMarker = getFSIMarkerID();
    iPoint = geometry->vertex[iMarker][iVertex];
    LoadDerivative = geometry->node[iPoint]->GetLoadDerivative();
    return LoadDerivative[0];
  }
  else {
    return 0.;
  }
}

double NativeSolidSolver::getLoadDerivativeY(unsigned short iVertex){
  if (config->GetKindProblem() == "ADJOINT")
  {
    unsigned short iMarker;
    unsigned long iPoint;
    double *LoadDerivative;
    iMarker = getFSIMarkerID();
    iPoint = geometry->vertex[iMarker][iVertex];
    LoadDerivative = geometry->node[iPoint]->GetLoadDerivative();
    return LoadDerivative[1];
  }
  else {
    return 0.;
  }
}

double NativeSolidSolver::getLoadDerivativeZ(unsigned short iVertex){
  if (config->GetKindProblem() == "ADJOINT")
  {
    unsigned short iMarker;
    unsigned long iPoint;
    double *LoadDerivative;
    iMarker = getFSIMarkerID();
    iPoint = geometry->vertex[iMarker][iVertex];
    LoadDerivative = geometry->node[iPoint]->GetLoadDerivative();
    return LoadDerivative[2];
  }
  else {
    return 0.;
  }
}

double NativeSolidSolver::getPlungeDampingDerivative()
{
  double dJdc = 0.;
  double hdot;
  unsigned short nInst = 2*getNumberHarmonics()+1;
  dJdc += integrator->GetSolver()->GetDampingDerivative(0);
  if (config->GetObjFunction() == "DISSIPATED_PLUNGE_POWER")
  {
    for (unsigned short iInst = 0; iInst < nInst; iInst++)
    {
      hdot = (integrator->GetSolver()->GetVel())[iInst];
      dJdc += hdot*hdot;
    }
  }
  return dJdc;
}

void NativeSolidSolver::computeInterfaceAdjointLoads(){

    double *Coord, Center[3], LoadDer[3];
    unsigned short iMarker, iVertex, nDim(3);
    unsigned long iPoint;
    CVector AdjLoads = integrator->GetSolver()->GetAdjLoads();

    /*--- Get the current center of rotation (can vary at each iteration) ---*/
    Center[0] = structure->GetCenterOfRotation_x();
    Center[1] = structure->GetCenterOfRotation_y();
    Center[2] = structure->GetCenterOfRotation_z();
    cout << "AdjLoads[0]: " << AdjLoads[0] << " AdjLoads[1]: " << AdjLoads[1] << endl;
    for(iMarker = 0; iMarker < geometry->GetnMarkers(); iMarker++){
      if (geometry->markersMoving[iMarker] == true){
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
          LoadDer[0] = 0.0;
          LoadDer[1] = 0.0;
          LoadDer[2] = 0.0;
          iPoint = geometry->vertex[iMarker][iVertex];
          Coord = geometry->node[iPoint]->GetCoord();
          LoadDer[0] += AdjLoads[1]*(Center[1]-Coord[1]); // Dependence of pitch eq. on X axis forces
          LoadDer[1] += AdjLoads[1]*(Center[0]-Coord[0]); // Dependence of pitch eq. on Y axis forces
          LoadDer[1] -= AdjLoads[0]; // Dependence of plunge eq. on Y axis forces
          geometry->node[iPoint]->SetLoadDerivative(LoadDer);
        }
      }
    }

}

void NativeSolidSolver::computeInterfaceAdjointLoads(unsigned int instance){

    double *Coord, Center[3], LoadDer[3];
    unsigned short iMarker, iVertex, nDim(3);
    unsigned int offset = instance+2*config->GetNumberHarmonics()+1; // Location in vector of pitch equation
    unsigned long iPoint;
    CVector AdjLoads = integrator->GetSolver()->GetAdjLoads();

    /*--- Get the current
     center of rotation (can vary at each iteration) ---*/
    Center[0] = structure->GetCenterOfRotation_x();
    Center[1] = structure->GetCenterOfRotation_y();
    Center[2] = structure->GetCenterOfRotation_z();
    cout << "Offset: " << offset << " AdjLoads[0]: " << AdjLoads[instance];
    if (structure->GetnDof()==2) cout << " AdjLoads[1]: " << AdjLoads[offset];
    cout << endl;
    for(iMarker = 0; iMarker < geometry->GetnMarkers(); iMarker++){
      if (geometry->markersMoving[iMarker] == true){
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
          LoadDer[0] = 0.0;
          LoadDer[1] = 0.0;
          LoadDer[2] = 0.0;
          iPoint = geometry->vertex[iMarker][iVertex];
          Coord = geometry->node[iPoint]->GetCoord();
          if (structure->GetnDof()==2){
            LoadDer[0] += AdjLoads[offset]*(Coord[1]-Center[1]); // -Dependence of pitch eq. on X axis forces
            LoadDer[1] -= AdjLoads[offset]*(Coord[0]-Center[0]); // -Dependence of pitch eq. on Y axis forces
          }
          LoadDer[1] -= AdjLoads[instance]; // -Dependence of plunge eq. on Y axis forces
          geometry->node[iPoint]->SetLoadDerivative(LoadDer);
        }
      }
    }

}
