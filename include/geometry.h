#pragma once

#include "config.h"
#include <vector>

class Point{

protected:
  double* Coord0;
  double* VarCoord; // Change in coordinates due to design variable application
  double* Coord;
  double* Coord_n;
  double* Vel;
  double* Vel_n;
  double* Force;
  double* DispAdj;
  double* LoadDerivative;

public:
  Point();
  ~Point();
  double* GetCoord0() const;
  double* GetCoord() const;
  double* GetVarCoord() const;
  double* GetCoord_n() const;
  double* GetVel() const;
  double* GetVel_n() const;
  double* GetForce() const;
  double* GetDisplacementAdjoint() const;
  double* GetLoadDerivative() const;
  void PrintCoord() const;
  void SetCoord0(double* newCoord);
  void SetVarCoord(double* newCoord);
  void SetCoord(double* newCoord);
  void SetCoord_n(double* newCoord);
  void SetVel(double* newVel);
  void SetVel_n(double* newVel_n);
  void SetForce(double* newForce);
  void SetDisplacementAdjoint(double* newDispAdj);
  void SetLoadDerivative(double* newLoadDerivative);
  void UpdateCoord();
  void UpdateVel();
};

/*class Vertex{

protected:
  double* Coord;

public:
  Vertex();
  ~Vertex();
  double* GetCoord() const;
  void SetCoord(double* newCoord);

}*/

/*class Element{

protected:
  Point** pointElem;
  unsigned short nPointElem;

public:
  Element();
  ~Element();
}*/

/*class Marker{

protected:
  std::string MarkerTag;
  unsigned long nElem;

public:
  Marker(std::string Tag, unsigned long nElem_value);
  ~Marker();
  Vertex** vertex;
  Element** element;

}*/

class Geometry{

protected:
  unsigned short nDim;
  unsigned long nElem, nPoint, nMarkers;
  unsigned long* nElemMarker;

public:
  Geometry(Config* config);
  ~Geometry();
  unsigned long GetnMarkers() const;
  bool GetMarkersMoving(unsigned long iMarker) const;
  unsigned long GetnPoint();
  void UpdateGeometry();
  unsigned long** vertex;  //vertex[iMarker][iPoint]
  unsigned long* nVertex;  //nVertex[iMarker]
  bool* markersMoving;
  Point** node;         //node[iPoint] returns a pointeur to a Point object so e.g. node[iPoint]->PrintCoord();
};

bool isInVec(std::vector<int> const& inputVector, int dummyInt);
void vecCopy(std::vector<int> const& sourceVector, unsigned long* destinationTab);
