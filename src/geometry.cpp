#include "../include/geometry.h"
#include "../include/config.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>

using namespace std;

Point::Point(){

  Coord0 = new double[3];
  Coord0[0] = 0.0;
  Coord0[1] = 0.0;
  Coord0[2] = 0.0;

  Coord = new double[3];
  Coord[0] = 0.0;
  Coord[1]= 0.0;
  Coord[2] = 0.0;

  Coord_n = new double[3];
  Coord_n[0] = 0.0;
  Coord_n[1] = 0.0;
  Coord_n[2] = 0.0;

  Vel = new double[3];
  Vel[0] = 0.0;
  Vel[1] = 0.0;
  Vel[2] = 0.0;

  Vel_n = new double[3];
  Vel_n[0] = 0.0;
  Vel_n[1] = 0.0;
  Vel_n[2] = 0.0;

  Force = new double[3];
  Force[0] = 0.0;
  Force[1] = 0.0;
  Force[2] = 0.0;

}

Point::~Point(){
  delete [] Coord0;
  Coord0 = NULL;

  delete [] Coord;
  Coord = NULL;

  delete [] Coord_n;
  Coord_n = NULL;

  delete [] Vel;
  Vel = NULL;

  delete [] Vel_n;
  Vel_n = NULL;

  delete [] Force;
  Force = NULL;
}

double* Point::GetCoord0() const{
  return Coord0;
}

double* Point::GetCoord() const{
  return Coord;
}

double* Point::GetCoord_n() const{
  return Coord_n;
}

double* Point::GetVel() const{
    return Vel;
}

double* Point::GetVel_n() const{
    return Vel_n;
}

double* Point::GetForce() const{
    return Force;
}

void Point::PrintCoord() const{
  cout << Coord[0] << " ; " << Coord[1] << " ; " << Coord[2] << endl;
}

void Point::SetCoord0(double* newCoord){
  Coord0[0] = newCoord[0];
  Coord0[1] = newCoord[1];
  Coord0[2] = newCoord[2];
}

void Point::SetCoord(double* newCoord){
  Coord[0] = newCoord[0];
  Coord[1] = newCoord[1];
  Coord[2] = newCoord[2];
}

void Point::SetCoord_n(double* newCoord){
  Coord_n[0] = newCoord[0];
  Coord_n[1] = newCoord[1];
  Coord_n[2] = newCoord[2];
}

void Point::SetVel(double *newVel){
    Vel[0] = newVel[0];
    Vel[1] = newVel[1];
    Vel[2] = newVel[2];
}

void Point::SetVel_n(double *newVel_n){
    Vel_n[0] = newVel_n[0];
    Vel_n[1] = newVel_n[1];
    Vel_n[2] = newVel_n[2];
}

void Point::SetForce(double* newForce){
    Force[0] = newForce[0];
    Force[1] = newForce[1];
    Force[2] = newForce[2];
}

void Point::UpdateCoord(){
  Coord_n[0] = Coord[0];
  Coord_n[1] = Coord[1];
  Coord_n[2] = Coord[2];
}

void Point::UpdateVel(){
    Vel_n[0] = Vel[0];
    Vel_n[1] = Vel[1];
    Vel_n[2] = Vel[2];
}

Geometry::Geometry(Config* config){

  string meshFileName, textLine, tampon;
  ifstream meshFile;
  string::size_type position;
  double Coord[3];
  double* TempCoord;
  unsigned long iMarker(0);
  int elemType(0), dummyInt(0);
  int iPoint;

  Coord[0] = 0.0;
  Coord[1] = 0.0;
  Coord[2] = 0.0;

  nDim = 0;
  nElem = 0;
  nPoint = 0;
  nMarkers = 0;
  //strcpy(cstr, config->GetMeshFile());

  meshFileName = config->GetMeshFile();

  /*--- Open the mesh file and check ---*/
  meshFile.open(meshFileName.c_str(), ios::in);
  if (meshFile.fail()){
    cout << "Unable to open the mesh file " << meshFileName << endl;
    exit(1);
  }

  cout << "Mesh file " << meshFileName << " is open." << endl;
  cout << "Constructing the mesh..." << endl;
  while(getline(meshFile, textLine)){

    position = textLine.find("NDIME=",0);
    if (position != string::npos){
      textLine.erase(0,6);
      nDim = atoi(textLine.c_str());
      cout << "Number of dimensions : " << nDim << endl;
    }

    position = textLine.find("NELEM=",0);
    if (position != string::npos){
      textLine.erase(0,6);
      nElem = atoi(textLine.c_str());
      cout << "Number of elements : " << nElem << endl;
      for(int iElem=0; iElem < nElem; iElem++){
        getline(meshFile, textLine);
      }
    }

    position = textLine.find("NPOIN=", 0);
    if (position != string::npos){
      textLine.erase(0,6);
      nPoint = atoi(textLine.c_str());
      cout << "Number of points : " << nPoint << endl;
      node = new Point*[nPoint];
      //cout << "JE VAIS REMPLIR" << endl;
      for(iPoint=0; iPoint < nPoint; iPoint++){
        getline(meshFile, textLine);
        //if(iPoint == 1) cout << textLine << endl;
        node[iPoint] = new Point();
        istringstream point_line(textLine);
        point_line >> Coord[0];
        //if(iPoint == 1) cout << Coord[0] << endl;
        point_line >> Coord[1];
        //if(iPoint == 1) cout << Coord[1] << endl;
        if(nDim == 3) point_line >> Coord[2];
        node[iPoint]->SetCoord0(Coord);
        node[iPoint]->SetCoord(Coord);
        node[iPoint]->SetCoord_n(Coord);
        TempCoord = node[iPoint]->GetCoord();
      //cout << iPoint << endl;
      }
    }

    position = textLine.find("NMARK=", 0);
    if (position != string::npos){
      textLine.erase(0,6);
      nMarkers = atoi(textLine.c_str());
      cout << "Number of markers : " << nMarkers << endl;
      vertex = new unsigned long*[nMarkers];
      nVertex = new unsigned long[nMarkers];
      markersMoving = new bool[nMarkers];
      nElemMarker = new unsigned long[nMarkers];
    }

    position = textLine.find("MARKER_TAG=",0);
    if (position != string::npos){
      textLine.erase(0,12);
      cout << textLine << endl;
      cout << "Reading elements for marker : " << textLine << endl;
      if (textLine.compare(config->GetMovingMarker()) == 0){
        cout << "Marker " << textLine << " is a moving marker." << endl;
        markersMoving[iMarker] = true;
      }
      else{
        markersMoving[iMarker] = false;
        cout << "Only marker " << config->GetMovingMarker() << "is a moving marker." << endl;
      } 
    }

    position = textLine.find("MARKER_ELEMS=",0);
    if (position != string::npos){
      textLine.erase(0,13);
      nElemMarker[iMarker] = atoi(textLine.c_str());
      cout << "Number of elements on the marker : " << nElemMarker[iMarker] << endl;
      vector<int> tempVertexMarker;
      for(int iElem=0; iElem < nElemMarker[iMarker]; iElem++){
        getline(meshFile,textLine);
        istringstream point_line(textLine);
        point_line >> elemType;
        if (elemType == 3){
          point_line >> dummyInt;
          if (!isInVec(tempVertexMarker, dummyInt)){
            tempVertexMarker.push_back(dummyInt);
          }
          point_line >> dummyInt;
          if (!isInVec(tempVertexMarker, dummyInt)){
            tempVertexMarker.push_back(dummyInt);
          }
        }
        else{
          cout << "Elem type " << elemType << " is not recognized !" << endl;
          exit(1);
        }
      }
      nVertex[iMarker] = tempVertexMarker.size();
      vertex[iMarker] = new unsigned long[nVertex[iMarker]];
      vecCopy(tempVertexMarker,vertex[iMarker]);
      iMarker++;
    }

  }

  meshFile.close();

  /*for (int iVertex = 0; iVertex < nVertex[0]; iVertex++){
    //cout << vertex[0][iVertex] << endl;
    iPoint = vertex[0][iVertex];
    node[iPoint]->PrintCoord();
  }*/



}

Geometry::~Geometry(){

  for(int iPoint=0; iPoint < nPoint; iPoint++){
    delete node[iPoint];
  }

  for(int iMarker=0; iMarker < nMarkers; iMarker++){
    delete [] vertex[iMarker];
  }

  delete [] vertex;

  delete [] node;

  delete [] nElemMarker;

}

unsigned long Geometry::GetnMarkers() const{
  return nMarkers;
}

bool Geometry::GetMarkersMoving(unsigned long iMarker) const{
  return markersMoving[iMarker];
}

void Geometry::UpdateGeometry(){

  unsigned long iPoint;

  for(iPoint=0; iPoint < nPoint; iPoint++){
    node[iPoint]->UpdateCoord();
    node[iPoint]->UpdateVel();
  }
}

bool isInVec(vector<int> const& inputVector, int dummyInt){
  int i = 0;
  while(i<inputVector.size() && inputVector[i]!=dummyInt) i++;
  return i<inputVector.size();
}

void vecCopy(vector<int> const& sourceVector, unsigned long* destinationTab){
  for(int i=0; i<sourceVector.size(); i++){
    destinationTab[i] = sourceVector[i];
  }
}
