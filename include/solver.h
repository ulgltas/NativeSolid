#pragma once

#include "MatVec.h"
#include "structure.h"
#include "config.h"
#include <iostream>

class Solver{

protected:
    CVector q;
    CVector qdot;
    CVector qddot;
    CVector q_n;
    CVector qdot_n;
    CVector qddot_n;
    CVector Loads;
    CVector Loads_n;
    CVector a;
    CVector a_n;
    bool linear;
    double L2norm;
    double dL2dwnorm;

public:
    Solver(unsigned int nDof, bool bool_linear);
    virtual ~Solver();
    virtual void Iterate(double& t0, double& tf, Structure* structure);
    virtual CVector & GetDisp();
    virtual CVector & GetVel();
    virtual CVector & GetAcc();
    virtual CVector & GetDisp_n();
    virtual CVector & GetVel_n();
    virtual CVector & GetAcc_n();
    virtual CVector & GetLoads();
    virtual CVector & GetAccVar();
    virtual CVector & GetAccVar_n();
    virtual void ResetSolution();
    virtual void SaveToThePast();
    virtual void SetInitialState(Config *config, Structure* structure);
    virtual void SetStateLoads(unsigned int iInstance, double load);
    virtual void SetStates(unsigned int iInstance, unsigned int dof,  double displacement);
    inline double GetL2Norm() {std::cout << L2norm << std::endl;return L2norm;}
    inline double GetdL2dwNorm() {return dL2dwnorm;}
    virtual void SetOmega(double val_omega) {};
    virtual double GetOmega() {return 0.;};

};

class AlphaGenSolver : public Solver {

protected:
    double beta;
    double gamma;
    double alpha_m;
    double alpha_f;
    double rho;
    double gammaPrime;
    double betaPrime;

public:
  AlphaGenSolver(unsigned int nDof, double val_rho, bool bool_linear);
  ~AlphaGenSolver();
  CVector & GetAccVar();
  CVector & GetAccVar_n();
  virtual void Iterate(double& t0, double& tf, Structure* structure);
  void ComputeRHS(Structure* structure, CVector &RHS);
  void ComputeResidual(Structure* structure, CVector &res);
  void ComputeTangentOperator(Structure* structure, CMatrix & St);
  void ResetSolution();
  void SaveToThePast();
  virtual void SetInitialState(Config *config, Structure* structure);

};

class RK4Solver : public Solver {

protected:
  unsigned int size;
  double lastTime;
  double currentTime;

public:
    RK4Solver(unsigned nDof, bool bool_linear);
    ~RK4Solver();
    virtual void Iterate(double &t0, double &tf, Structure* structure);
    void EvaluateStateDerivative(double tCurrent, CVector& state, CVector& stateDerivative, Structure* structure);
    void interpLoads(double& tCurrent, CVector& val_loads);
    virtual void SetInitialState(Config* config, Structure* structure);
    CVector SetState();
    CVector SetState_n();

};

class StaticSolver : public Solver {

protected:
    unsigned int _nDof;
    CMatrix KK;

public:
    StaticSolver(unsigned nDof, bool bool_linear);
    ~StaticSolver();

    virtual void Iterate(double &t0, double &tf, Structure* structure);
    virtual void SetInitialState(Config* config, Structure* structure);

};

class HarmonicSolver : public Solver {
protected:
    unsigned int _nHarmonic; // Number of Harmonics
    unsigned int _nOmega; // Number of Harmonics*2 + 1
    unsigned int _nDof; // So far one DoF
    double omega;
    CMatrix d; // Harmonic balance time derivative
    CMatrix d2; // Harmonic balance second time derivative
    CMatrix AA; // Harmonic balance physics matrix
    CMatrix  E; // Harmonic balance DFT matrix
    CMatrix Em1; // Harmonic balance IFT matrix
    CVector stateLoads; // Load at each time interval

public:
    HarmonicSolver(unsigned nDof, unsigned int nHarmonic, bool bool_linear);
    ~HarmonicSolver();

    virtual void Iterate(double& t0, double& tf, Structure* structure);
    virtual void SetInitialState(Config* config, Structure* structure);
    virtual void SetHBMatrices();
    virtual void SetStateLoads(unsigned int iInstance, double load);
    virtual void SetStates(unsigned int iInstance, unsigned int dof,  double displacement);
    virtual void SetOmega(double val_omega) {omega = val_omega; SetHBMatrices();};
    virtual double GetOmega() {return omega;};
    /*inline double GetL2Norm() {std::cout << L2norm << std::endl; return L2norm;}
    inline double GetdL2dwNorm() {return dL2dwnorm;}*/
};
