#include "../include/MatVec.h"
#include <iostream>
#include <cmath>

using namespace std;

/*--- Definition of the CVector class ---*/

CVector::CVector(void):nElm(0){

  vec_val = NULL;

}

CVector::CVector(const unsigned long & size, const double & val){

  vec_val = NULL;
  if(size <= 0){
    nElm = 0;
    cerr << "CVector:CVector(const unsigned long &, const double): " << "invalid input : size = " << size << endl;
    throw(-1);
  }
  else{
    nElm = size;
    vec_val = new double[nElm];
    for(unsigned int i=0; i < nElm; i++)
      vec_val[i] = val;
  }
}

CVector::CVector(const CVector & u){

  nElm = u.nElm;
  vec_val = NULL;
  vec_val = new double[nElm];
  for(unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];
}

CVector::~CVector(){

  if (vec_val != NULL) delete [] vec_val;
}

void CVector::Initialize(const unsigned long &size, const double & val){

  if(vec_val == NULL){
    if(size <= 0){
      nElm = 0;
      cerr << "CVector::Initialize(const unsigned long &, const double &): " << "invalid number of element" << nElm << endl;
      throw(-1);
    }
    else{
      nElm = size;
      vec_val = new double[nElm];
      for(unsigned int ii=0; ii<nElm; ii++){
        vec_val[ii] = val;
      }
    }
  }
  else{
    std::cerr << "CVector::Initialize(const unsigned long &, const double &): " << "vector already initialize, resizing is not allowed." << std::endl;
    throw(-1);
  }
}

unsigned long CVector::GetSize() const{

  if(nElm <= 0){
    cerr << "CVector::GetSize() const: " << "invalid number of element" << nElm << endl;
    throw(-1);
  }

  return nElm;
}

double* CVector::GetVec() const{
  return vec_val;
}

void CVector::print() const{
  cout << "**********" << endl;
  for(unsigned int i=0; i<nElm; i++)
    cout << vec_val[i] << endl;
  cout << "**********" << endl;
}

void CVector::SetAllValues(const double & val){

  for(unsigned int ii=0; ii<nElm; ii++){
    vec_val[ii] = val;
  }
}

CVector & CVector::operator=(const CVector & u){
  if(this == &u) return *this;

  if(nElm != u.nElm){
     cerr << "Sizes do not match for operator = : size1 = " << nElm << " & size2 = " << u.nElm << endl;
     throw(-1);
  }
  else{
    for(unsigned int i=0; i<nElm; i++)
      vec_val[i] = u.vec_val[i];
  }

  return *this;
}

CVector & CVector::operator+=(const CVector & u){
   if(nElm != u.nElm){
    cerr << "CVector::operator+=(const CVector &) const: " << "sizes do not match" << endl;
    throw(-1);
   }

    for(unsigned int i=0; i < nElm; i++)
        vec_val[i] += u.vec_val[i];

    return *this;
 }

CVector & CVector::operator-=(const CVector & u){
   if(nElm != u.nElm){
    cerr << "CVector::operator-=const CVector &) const: " << "sizes do not match" << endl;
    throw(-1);
   }

   for(unsigned int i=0; i<nElm; i++)
    vec_val[i] -= u.vec_val[i];

   return *this;
 }

CVector & CVector::operator*=(const double & val){
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] *= val;

  return *this;
}

CVector & CVector::operator/=(const double & val){
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] /= val;

  return *this;
}

double & CVector::operator[](const unsigned long & i) const{
  return vec_val[i];
}

double CVector::dotProd(const CVector & v) const{
  if(nElm != v.nElm){
    cerr << "CVector::dotProd(const CVector &, const CVector &): " << "sizes do not math" << endl;
    throw(-1);
  }

  double dotProd(0.0);
  for(unsigned int i=0; i<nElm; i++)
    dotProd += vec_val[i]*v.vec_val[i];

  return dotProd;
}

double CVector::norm() const{
  double norm(0), tempVal(0);
  tempVal = (*this).dotProd(*this);

  if(tempVal < 0){
    cerr << "CVector::norm(): " << "computed dot product is negative: " << tempVal << endl;
    throw(-1);
  }

  norm = sqrt(tempVal);

  return norm;
}

void CVector::Reset(){
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] = 0;
}

/*--- Definition of the CMatrix class ---*/

CMatrix::CMatrix(void):nEq(0), nVar(0){}

CMatrix::CMatrix(const unsigned long & val_nEq, const unsigned long & val_nVar, const double & val){
  nEq = val_nEq;
  nVar = val_nVar;
  unsigned long nElm = nEq*nVar;

  mat_val.Initialize(nElm, val);
}

CMatrix::CMatrix(const CMatrix & A){

  nEq = A.nEq;
  nVar = A.nVar;
  unsigned long nElm = nEq*nVar;

  mat_val.Initialize(nElm, 0.0);
  mat_val = A.mat_val;
}

CMatrix::~CMatrix(){}

void CMatrix::Initialize(const unsigned long & val_nEq, const unsigned long & val_nVar, const double & val){

  nEq = val_nEq;
  nVar = val_nVar;
  unsigned long nElm = nEq*nVar;

  mat_val.Initialize(nElm, val);
}

void CMatrix::print() const{
  mat_val.print();
}

CMatrix & CMatrix::operator=(const CMatrix & a){
  if(this == &a) return *this;

  if(nEq != a.nEq || nVar != a.nVar){
    cerr << "Sizes do not match for operator = : nEq1 = " << nEq << " & nEq2 = " << a.nEq << " ; nVar1 = " << nVar << " & nVar2 = " << a.nVar << endl;
    throw(-1);
  }
  else{
    mat_val = a.mat_val;
  }

  return *this;
}

CMatrix & CMatrix::operator+=(const CMatrix & a){

  if(nEq != a.nEq || nVar != a.nVar){
    cerr << "CMatrix::operator+=(const CMatrix &): " << "sizes do not match" << endl;
    throw(-1);
  }
  else{
    mat_val += a.mat_val;
  }

  return *this;
}

CMatrix & CMatrix::operator-=(const CMatrix & a){

  if(nEq != a.nEq || nVar != a.nVar){
    cerr << "CMatrix::operator+=(const CMatrix &): " << "sizes do not match" << endl;
    throw(-1);
  }
  else{
    mat_val -= a.mat_val;
  }

  return *this;
}

CMatrix & CMatrix::operator*=(const double & val){

  mat_val *= val;

  return *this;

}

CMatrix & CMatrix::operator/=(const double & val){

  mat_val /= val;

  return *this;

}

unsigned long CMatrix::GetnEq() const{
  return nEq;
}

unsigned long CMatrix::GetnVar() const{
  return nVar;
}

double* CMatrix::GetMat() const{
  return mat_val.GetVec();
}

CVector CMatrix::GetCVec() const{
  return mat_val;
}

void CMatrix::SetElm(const int &i, const int &j, double val){
  mat_val[(j-1)*nEq+(i-1)] = val;
}

double CMatrix::GetElm(const int &i, const int &j) const{
  return mat_val[(j-1)*nEq+(i-1)];
}

double CMatrix::DiagProduct() const{
  double prod(1);
  if (nEq == nVar){
    for(unsigned int i=1; i<=nEq; i++){
      prod *= GetElm(i,i);
    }
  }

  return prod;
}

double CMatrix::ComputeDet() const{
  double Det(1);
  int NbChange(0);
  double PLU[nEq*nVar];
  int IPIV[nEq];
  int INFO(-1);

  *PLU = *(mat_val.GetVec());
  INFO = LAPACKE_dgetrf(LAPACK_COL_MAJOR,nEq,nVar,PLU,nEq,IPIV);

  for(int i=0; i<nEq-1; i++){
    if(IPIV[i] != i+1) NbChange++;
  }

  return 0.0;
}

void CMatrix::Reset(){
  mat_val.Reset();
}

CVector MatVecProd(const CMatrix &A, const CVector &b){

  CVector x(A.GetnVar());
  cblas_dgemv(CblasColMajor,CblasNoTrans,A.GetnEq(),A.GetnVar(),1,A.GetMat(),A.GetnEq(),b.GetVec(),1,0,x.GetVec(),1);
  return x;

}

CMatrix MatMatProd(const CMatrix &A, const CMatrix &B){

  CMatrix X(A.GetnEq(), B.GetnVar(), 0);
  if (A.GetnVar() != B.GetnEq())
  {
    cerr << "CMatrix MatMatProd(const CMatrix &A, const CMatrix &B): " << "invalid inputs : size of A = " << A.GetnEq() << "x" << A.GetnVar() << "size of B = " << B.GetnEq() << "x" << B.GetnVar()<< endl;
    throw(-1);
  }
  
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,A.GetnEq(),B.GetnVar(),A.GetnVar(),1.0,A.GetMat(),A.GetnEq(),B.GetMat(),B.GetnEq(),0.0,X.GetMat(),X.GetnEq());
  return X;

}

CVector ScalVecProd(const double & scal, const CVector &b){

  CVector x(b.GetSize());
  cblas_daxpy(b.GetSize(),scal,b.GetVec(),1,x.GetVec(),1);
  return x;

}

CMatrix ScalMatProd(const double & scal, const CMatrix &A){

  CMatrix copy(A);
  copy *= scal;
  return copy;
}

/* INDEPENDENT FUNCTIONS */
int SolveSys(const CMatrix &A, CVector &b){
  const int N = A.GetnEq();
  const int LDA = A.GetnVar();
  CMatrix Mat(A);
  int NRHS(1);
  int *IPIV;
  int LDB = b.GetSize();
  int INFO(-1);

  IPIV = new int[N];

  INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR,N,NRHS,Mat.GetMat(),LDA,IPIV,b.GetVec(),LDB);

  return INFO;
}

int SolveSys(const CMatrix &A, CMatrix &B){
  const int N = A.GetnEq();
  const int LDA = A.GetnVar();
  CMatrix Mat(A);
  int NRHS = B.GetnEq();
  int *IPIV;
  int LDB = B.GetnVar();
  int INFO(-1);

  IPIV = new int[N];

  INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR,N,NRHS,Mat.GetMat(),LDA,IPIV,B.GetMat(),LDB);

  return INFO;
}

void MatrixToVec(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec){
    if(order == COL_MAJ){
      for(int j=0; j<sizeVec; j++)
        vecteur[j] = matrix[j%Nrow][j/Nrow];
    }
    else if(order == ROW_MAJ){
      for(int j=0; j<sizeVec; j++)
        vecteur[j] = matrix[j/Ncol][j%Ncol];
    }
    else{
      cerr << "Wrong storage order" << endl;
      throw(-1);
    }
}

void VecToMatrix(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec){
    if(order == COL_MAJ){
      for(int i=0; i<Nrow;i++)
        for(int j=0; j<Ncol; j++)
          matrix[i][j] = vecteur[j*Nrow+i];
    }
    else if(order == ROW_MAJ){
      for(int i=0; i<Nrow;i++)
        for(int j=0; j<Ncol; j++)
          matrix[i][j] = vecteur[i*Ncol+j];
    }
    else{
      cerr << "Wrong storage order" << endl;
      throw(-1);
    }
}

/* OPERATORS OVERLOADS*/
CVector operator+(CVector & vecA, CVector & vecB){

  CVector copy(vecA);
  copy += vecB;
  return copy;
}

CVector operator+(const CVector & vecA, const CVector & vecB){

  CVector copy(vecA);
  copy += vecB;
  return copy;
}

CVector operator+(CVector & vecA, const CVector & vecB){

  CVector copy(vecA);
  copy += vecB;
  return copy;
}

CVector operator+(const CVector & vecA, CVector & vecB){

  CVector copy(vecA);
  copy += vecB;
  return copy;
}



CVector operator-(CVector & vecA, CVector & vecB){

  CVector copy(vecA);
  copy -= vecB;
  return copy;
}

CVector operator-(const CVector & vecA, CVector & vecB){

  CVector copy(vecA);
  copy -= vecB;
  return copy;
}

CVector operator-(CVector & vecA, const CVector & vecB){

  CVector copy(vecA);
  copy -= vecB;
  return copy;
}

CVector operator-(const CVector & vecA, const CVector & vecB){

  CVector copy(vecA);
  copy -= vecB;
  return copy;
}



CVector operator*(CVector & vecA, double & val_mult){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(CVector & vecA, const double & val_mult){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(double & val_mult, CVector & vecA){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(const double & val_mult, CVector & vecA){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(const CVector & vecA, double & val_mult){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(const CVector & vecA, const double & val_mult){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(double & val_mult, const CVector & vecA){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}

CVector operator*(const double & val_mult, const CVector & vecA){

  CVector copy(vecA);
  copy *= val_mult;
  return copy;
}



CVector operator/(CVector & vecA, double & val_mult){

  CVector copy(vecA);
  copy /= val_mult;
  return copy;
}

CVector operator/(CVector & vecA, const double & val_mult){

  CVector copy(vecA);
  copy /= val_mult;
  return copy;
}

CVector operator/(const CVector & vecA, double & val_mult){

  CVector copy(vecA);
  copy /= val_mult;
  return copy;
}

CVector operator/(const CVector & vecA, const double & val_mult){

  CVector copy(vecA);
  copy /= val_mult;
  return copy;
}


CMatrix operator+(CMatrix & matA, CMatrix & matB){

  CMatrix copy(matA);
  copy += matB;
  return copy;
}

CMatrix operator-(CMatrix & matA, CMatrix & matB){

  CMatrix copy(matA);
  copy -= matB;
  return copy;
}

CMatrix operator*(CMatrix matA, double & val){

  CMatrix copy(matA);
  copy *= val;
  return copy;
}

CMatrix operator*(double & val, CMatrix matA){

  CMatrix copy(matA);
  copy *= val;
  return copy;
}
