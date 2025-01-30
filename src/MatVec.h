#pragma once

#include <cblas.h>
#include <lapacke.h>

#define ROW_MAJ 0
#define COL_MAJ 1

class CVector
{
    unsigned long nElm;
    double *vec_val;

public:
    CVector();
    CVector(const unsigned long &size, const double &val = 0.0);
    CVector(const CVector &u);
    ~CVector();
    unsigned long GetSize() const;
    double *GetVec() const;
    void print() const;
    void Initialize(const unsigned long &size, const double &val = 0.0);
    void SetAllValues(const double &val);
    CVector &operator=(const CVector &u);
    CVector &operator+=(const CVector &u);
    CVector &operator-=(const CVector &u);
    CVector &operator*=(const double &val);
    CVector &operator/=(const double &val);
    double &operator[](const unsigned long &i) const;
    double norm() const;
    double dotProd(const CVector &v) const;
    void Reset();
};

class CMatrix
{
    unsigned long nEq;
    unsigned long nVar;
    CVector mat_val;

public:
    CMatrix();
    CMatrix(const unsigned long &val_nEq, const unsigned long &val_nVar, const double &val = 0.0);
    CMatrix(const CMatrix &A);
    ~CMatrix();
    void Initialize(const unsigned long &val_nEq, const unsigned long &val_nVar, const double &val = 0.0);
    void print() const;
    CMatrix &operator=(const CMatrix &a);
    CMatrix &operator+=(const CMatrix &a);
    CMatrix &operator-=(const CMatrix &a);
    CMatrix &operator*=(const double &val);
    CMatrix &operator/=(const double &val);
    unsigned long GetnEq() const;
    unsigned long GetnVar() const;
    double *GetMat() const;
    CVector GetCVec() const;
    void SetElm(const int &i, const int &j, double val);
    double GetElm(const int &i, const int &j) const;
    double DiagProduct() const;
    double ComputeDet() const;
    void Reset();
};

CVector MatVecProd(const CMatrix &A, const CVector &b);
CVector ScalVecProd(const double &scal, const CVector &b);
CMatrix ScalMatProd(const double &scal, const CMatrix &A);
CMatrix MatMatProd(const CMatrix &A, const CMatrix &B);
int SolveSys(const CMatrix &A, CVector &b);
int SolveSys(const CMatrix &A, CMatrix &B);

void MatrixToVec(int order, double **matrix, double *vecteur, int Nrow, int Ncol, int sizeVec);
void VecToMatrix(int order, double **matrix, double *vecteur, int Nrow, int Ncol, int sizeVec);

// OPERATOR +
CVector operator+(CVector &vecA, CVector &vecB);
CVector operator+(const CVector &vecA, const CVector &vecB);
CVector operator+(CVector &vecA, const CVector &vecB);
CVector operator+(const CVector &vecA, CVector &vecB);
CMatrix operator+(const CMatrix &matA, const CMatrix &matB);

// OPERATOR -
CVector operator-(CVector &vecA, CVector &vecB);
CVector operator-(const CVector &vecA, CVector &vecB);
CVector operator-(CVector &vecA, const CVector &vecB);
CVector operator-(const CVector &vecA, const CVector &vecB);
CMatrix operator-(const CMatrix &matA, const CMatrix &matB);

// OPERATOR *
CVector operator*(CVector &vecA, double &val_mult);
CVector operator*(CVector &vecA, const double &val_mult);
CVector operator*(const CVector &vecA, double &val_mult);
CVector operator*(const CVector &vecA, const double &val_mult);
CVector operator*(double &val_mult, CVector &vecA);
CVector operator*(const double &val_mult, CVector &vecA);
CVector operator*(double &val_mult, const CVector &vecA);
CVector operator*(const double &val_mult, const CVector &vecA);
CMatrix operator*(CMatrix matA, double &val);
CMatrix operator*(double &val, CMatrix matA);

// OPERATOR /
CVector operator/(CVector &vecA, double &val_mult);
CVector operator/(CVector &vecA, const double &val_mult);
CVector operator/(const CVector &vecA, double &val_mult);
CVector operator/(const CVector &vecA, const double &val_mult);
