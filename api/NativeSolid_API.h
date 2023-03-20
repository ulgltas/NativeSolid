#pragma once

#include <string>
#include "../include/config.h"
#include "../include/structure.h"
#include "../include/integration.h"
#include "../include/output.h"
#include "../include/MatVec.h"
#include "../include/geometry.h"
#include <fstream>
#include <iostream>
#include <sstream>


class NativeSolidSolver{

protected:
    std::string confFile;
    Config* config;
    Geometry* geometry;
    Structure* structure;
    Integration* integrator;
    Output* output;
    std::ofstream historyFile;
    std::ofstream historyFile2;
    std::ofstream restartFile;
    CVector q_uM1; //The displacement at the previous FSI iteration
    double omega;
    unsigned long nSolidInterfaceVertex;
    double varCoordNorm;


public:
    /* NEW generation */
    NativeSolidSolver(std::string str, bool FSIComp);
    ~NativeSolidSolver();
    void exit();
    //double getVarCoordNorm() const;
    void preprocessIteration(unsigned long ExtIter);
    void timeIteration(double t0, double tf);
    //void mapRigidBodyMotion(bool predicition, bool initialize);
    void computeInterfacePosVel(bool initialize);
    void computeInterfacePosVel(bool initialize, unsigned int instance);
    void computeInterfaceAdjointLoads();
    void computeInterfaceAdjointLoads(unsigned int instance);
    void setInitialDisplacements();
    void staticComputation();
    void harmonicComputation();
    void writeSolution(double currentTime, double lastTime, double currentFSIIter, unsigned long ExtIter, unsigned long NbExtIter);
    void writeSolution(double time, int FSIter);
    void writeAdjointSolution(double time, int FSIter);
    void saveSolution();
    void updateSolution();
    //void updateGeometry();
    //void displacementPredictor();
    unsigned short getFSIMarkerID();
    unsigned long getNumberOfSolidInterfaceNodes(unsigned short iMarker);
    unsigned int getInterfaceNodeGlobalIndex(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosZ(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosX0(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosY0(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosZ0(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeDispX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeDispY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeDispZ(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelZ(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelXNm1(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelYNm1(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelZNm1(unsigned short iMarker, unsigned short iVertex);
    double getRotationCenterPosX();
    double getRotationCenterPosY();
    double getRotationCenterPosZ();
    void setGeneralisedForce(unsigned int instance = 0);
    void setGeneralisedForce(double Fx, double Fy);
    void setGeneralisedMoment(unsigned int instance = 0);
    void setGeneralisedMoment(double M);
    void applyload(unsigned short iVertex, unsigned int iInst, double Fx, double Fy, double Fz);
    unsigned int getNumberHarmonics();
    inline double getDeltaOmega() {return integrator->GetSolver()->GetDeltaOmega();};
    inline void setOmega(double val_omega) {integrator->GetSolver()->SetOmega(val_omega);}
    void applyDisplacementAdjoint(unsigned short iVertex, unsigned int instance, double Dx, double Dy, double Dz);
    double getLoadDerivativeX(unsigned short iVertex);
    double getLoadDerivativeY(unsigned short iVertex);
    double getLoadDerivativeZ(unsigned short iVertex);
    void setTotalAdjointDisplacement(unsigned int instance);
    void setFrequencyDerivative(double dJdw);
    double getObjectiveFunction();
};
