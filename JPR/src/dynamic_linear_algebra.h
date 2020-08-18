#ifndef dynamic_linear_algebra_h_
#define dynamic_linear_algebra_h_

#include <math.h>


template <typename FloatType>
class DynamicVector {
 public:
  //  ---- constructors, copy constructor
  inline DynamicVector<FloatType>(int dim = 0);
  inline DynamicVector<FloatType>(const DynamicVector<FloatType> &init);

  // ---- destructor
  inline ~DynamicVector();

  inline int GetDim() const;
  // alias to getDim()
  inline int Size() const;
  inline void SetDim(const int dim);
  // set all entries to (FloatType)0
  inline void SetZero();

  inline FloatType& operator[](const int &index);
  inline const FloatType& operator[](const int &index) const;

  //  ----  operators +, -, *, /, etc.
  inline DynamicVector<FloatType> operator+(const DynamicVector<FloatType> &op) const;
  inline DynamicVector<FloatType> operator-(const DynamicVector<FloatType> &op) const;
  inline DynamicVector<FloatType> operator-() const;
  inline FloatType operator*(const DynamicVector<FloatType> op) const;
  inline DynamicVector<FloatType> operator*(const FloatType &s) const;
  inline DynamicVector<FloatType> operator/(const FloatType &s) const;
  inline FloatType ComponentSum() const;
  inline FloatType Average() const;

  //  ---- operators +=, -=, *=, /=
  inline DynamicVector<FloatType> operator+=(const DynamicVector<FloatType> &op);
  inline DynamicVector<FloatType> operator-=(const DynamicVector<FloatType> &op);
  inline DynamicVector<FloatType> operator*=(const FloatType &s);
  inline DynamicVector<FloatType> operator/=(const FloatType &s);

  //  ---- operators =, ==, !=
  inline void operator=(const DynamicVector<FloatType> &op);
  inline bool operator==(const DynamicVector<FloatType> op) const;
  inline bool operator!=(const DynamicVector<FloatType> op) const;

  inline FloatType MaxComponent() const;
  inline FloatType MinComponent() const;
  inline DynamicVector<FloatType> AbsoluteValues() const;
 private:
  int dim;
  FloatType *v;
};

template <typename FloatType>
class DynamicMatrix {
 public:
  // ---- constructors, copy constructor
  inline DynamicMatrix( int c = 0, int r = 0, bool initAsIdentity = true);
  inline DynamicMatrix(const DynamicMatrix<FloatType> &init);

  // ---- destructor
  inline ~DynamicMatrix();

  inline void SetDimension(int num_columns_ = 0, int num_rows_ = 0, bool initAsIdentity = true);
  inline void SetZero();
  inline int GetRowsDim() const;
  inline int GetColsDim() const;
  inline int GetRows() const;
  inline int GetColumns() const;

  inline DynamicVector<FloatType>& operator[](const int &index);
  inline const DynamicVector<FloatType>& operator[](const int &index) const;

  //  ----  operators +, -, *, /
  inline DynamicMatrix<FloatType> operator+(const DynamicMatrix<FloatType> &op) const;
  inline DynamicMatrix<FloatType> operator-(const DynamicMatrix<FloatType> &op) const;
  inline DynamicMatrix<FloatType> operator-() const;
  inline DynamicMatrix<FloatType> operator*(const FloatType &s) const;
  inline DynamicMatrix<FloatType> operator*(const DynamicMatrix<FloatType> &op) const;
  inline DynamicMatrix<FloatType> operator/(const FloatType &s) const;

  //  ---- operators +=, -=, *=, /=
  inline DynamicMatrix<FloatType> operator+=(const DynamicMatrix<FloatType> &op);
  inline DynamicMatrix<FloatType> operator-=(const DynamicMatrix<FloatType> &op);
  inline DynamicMatrix<FloatType> operator*=(const FloatType &op);
  inline DynamicMatrix<FloatType> operator*=(const DynamicMatrix<FloatType> &op);
  inline DynamicMatrix<FloatType> operator/=(const FloatType &op);

  //  ---- operators =, ==, !=
  inline void operator=(const DynamicMatrix<FloatType> &op);
  inline bool operator==(const DynamicMatrix<FloatType> &op) const;
  inline bool operator!=(const DynamicMatrix<FloatType> &op) const;

  // ---- multiplying with vectors
  inline DynamicVector<FloatType> operator*(const DynamicVector<FloatType> &v) const;

  // --- Norm---
  inline FloatType FrobeniusNorm() const;

  inline DynamicMatrix<FloatType> Transpose() const;

 private:
  int num_rows_;
  int num_columns_;
  DynamicVector<FloatType> *theColumns;
};


// ---------------------------------------------------------------------------
//
//                                   instances
//
// ---------------------------------------------------------------------------


typedef DynamicVector<float> DVectorF;
typedef DynamicVector<double> DVectorD;
typedef DynamicVector<int> DVectorI;

typedef DynamicMatrix<float> DMatrixF;
typedef DynamicMatrix<double> DMatrixD;
typedef DynamicMatrix<int> DMatrixI;

#endif
