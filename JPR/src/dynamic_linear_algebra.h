#ifndef dynamic_linear_algebra_h_
#define dynamic_linear_algebra_h_

#include <math.h>


template <typename FloatType>
class DynamicVector {
 public:
  //  ---- constructors, copy constructor
  inline DynamicVector<FloatType>(int dim = 0) {  
      this->dim = dim;
      v = new FloatType[dim];
  }
  inline DynamicVector<FloatType>(const DynamicVector<FloatType> &init) {
      this->dim = init.dim;
      v = new FloatType[dim];
      for (int i = 0; i < dim; i++) {
          v[i] = init.v[i];
      }
  }

  // ---- destructor
  inline ~DynamicVector() {
      delete[] v;
  }

  inline int GetDim() const {
      return dim;
  }
  
  // alias to getDim()
  inline int Size() const {
      return dim;
  }
  inline void SetDim(const int dim) {
      this->dim = dim;
      delete[] v;
      v = new FloatType[dim];
  }
  // set all entries to (FloatType)0
  inline void SetZero() {
      for (int i = 0; i < dim; i++) {
          v[i] = (FloatType)0;
      }
  }

  inline FloatType& operator[](const int &index) {
      return v[index];
  }
  inline const FloatType& operator[](const int &index) const {
      return v[index];
  }

  //  ----  operators +, -, *, /, etc.
  inline DynamicVector<FloatType> operator+(const DynamicVector<FloatType> &op) const {
      DynamicVector<FloatType> result(dim);
      for (int i = 0; i < dim; i++) {
          result[i] = v[i] + op.v[i];
      }
      return result;
  }
  inline DynamicVector<FloatType> operator-(const DynamicVector<FloatType> &op) const {
      DynamicVector<FloatType> result(dim);
      for (int i = 0; i < dim; i++) {
          result[i] = v[i] - op.v[i];
      }
      return result;
  }
  inline DynamicVector<FloatType> operator-() const {
      DynamicVector<FloatType> result(dim);
      for (int i = 0; i < dim; i++) {
          result[i] = -v[i];
      }
      return result;
  }
  inline FloatType operator*(const DynamicVector<FloatType> op) const {
      FloatType result = 0.0;
      for (int i = 0; i < dim; i++) {
          result += v[i] * op.v[i];
      }
      return result;
  }
  inline DynamicVector<FloatType> operator*(const FloatType &s) const {
      DynamicVector<FloatType> result(dim);
      for (int i = 0; i < dim; i++) {
          result[i] = s * v[i];
      }
      return result;
  }
  inline DynamicVector<FloatType> operator/(const FloatType &s) const {
      DynamicVector<FloatType> result(dim);
      for (int i = 0; i < dim; i++) {
          result[i] = v[i] / s;
      }
      return result;
  }
  inline FloatType ComponentSum() const {
      FloatType result = 0;
      for (int i = 0; i < dim; i++) {
          result += v[i];
      }
      return result;
  }
  inline FloatType Average() const {
      return ComponentSum() / dim;  
  }

  //  ---- operators +=, -=, *=, /=
  inline DynamicVector<FloatType> operator+=(const DynamicVector<FloatType> &op) {
      for (int i = 0; i < dim; i++) {
          v[i] += op.v[i];
      }
      return *this;
  }
  inline DynamicVector<FloatType> operator-=(const DynamicVector<FloatType> &op) {
      for (int i = 0; i < dim; i++) {
          v[i] -= op.v[i];
      }
      return *this;  
  }
  inline DynamicVector<FloatType> operator*=(const FloatType &s) {
      for (int i = 0; i < dim; i++) {
          v[i] *= s;
      }
      return *this;
  }
  inline DynamicVector<FloatType> operator/=(const FloatType &s) {
      for (int i = 0; i < dim; i++) {
          v[i] /= s;
      }
      return *this;
  }

  //  ---- operators =, ==, !=
  inline void operator=(const DynamicVector<FloatType> &op) {
      if (dim != op.dim) {
          delete[] v;
          v = new FloatType[op.dim];
      }
      this->dim = op.dim;
      for (int i = 0; i < op.dim; i++) {
          v[i] = op.v[i];
      }
  }
  inline bool operator==(const DynamicVector<FloatType> op) const {
      for (int i = 0; i < dim; i++) {
          if (v[i] != op.v[i]) return false;
      }
      return true;
  }
  inline bool operator!=(const DynamicVector<FloatType> op) const {
      for (int i = 0; i < dim; i++) {
          if (v[i] != op.v[i]) return true;
      }
      return false;
  }

  inline FloatType MaxComponent() const {
      FloatType result = v[0];
      for (int i = 1; i < dim; i++) {
          if (v[i] > result) {
              result = v[i];
          }
      }
      return result;
  }
  inline FloatType MinComponent() const {
      FloatType result = v[0];
      for (int i = 1; i < dim; i++) {
          if (v[i] < result) {
              result = v[i];
          }
      }
      return result;
  }
  inline DynamicVector<FloatType> AbsoluteValues() const {
      DynamicVector<FloatType> result(dim);
      for (int i = 0; i < dim; i++) {
          result[i] = fabs(v[i]);
      }
      return result;
  }
 private:
  int dim;
  FloatType *v;
};



typedef DynamicVector<float> DVectorF;
typedef DynamicVector<double> DVectorD;
typedef DynamicVector<int> DVectorI;


// ---------------------------------------------------------------------------
//
//                                   instances
//
// --------------------------------------------------------------------------


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

template <typename FloatType>
inline DynamicMatrix<FloatType>::DynamicMatrix(int c, int r, bool initAsIdentity) {
  num_columns_ = c;
  num_rows_ = r;
  theColumns = new DynamicVector<FloatType>[num_columns_];
  for (int n = 0; n < num_columns_; n++)
    theColumns[n] = DynamicVector<FloatType>(num_rows_);
  for (int j = 0; j < num_columns_; j++) {
    for (int i = 0; i < num_rows_; i++) {
      if (i == j && initAsIdentity) {
        theColumns[j][i] = 1;
      }
      else {
        theColumns[j][i] = 0;
      }
    }
  }
}

template <typename FloatType>
inline void DynamicMatrix<FloatType>::SetDimension(int num_columns_, int num_rows_, bool initAsIdentity) {
  delete[] theColumns;
  this->num_rows_ = num_rows_;
  this->num_columns_ = num_columns_;
  theColumns = new DynamicVector<FloatType>[num_columns_];
  for (int n = 0; n < num_columns_; n++)
    theColumns[n] = DynamicVector<FloatType>(num_rows_);

  for (int j = 0; j < num_columns_; j++) {
    for (int i = 0; i < num_rows_; i++) {
      if (i == j && initAsIdentity) {
        theColumns[j][i] = 1;
      }
      else {
        theColumns[j][i] = 0;
      }
    }
  }
}

template <typename FloatType>
inline int DynamicMatrix<FloatType>::GetRowsDim() const {
  return num_rows_;
}

template <typename FloatType>
inline int DynamicMatrix<FloatType>::GetColsDim() const {
  return num_columns_;
}

template <typename FloatType>
inline void DynamicMatrix<FloatType>::SetZero() {
  for (int j = 0; j < num_columns_; j++)
    for (int i = 0; i < num_rows_; i++)
      theColumns[j][i] = 0.0;
}

template <typename FloatType>
inline int DynamicMatrix<FloatType>::GetRows() const {
  return num_rows_;
}

template <typename FloatType>
inline int DynamicMatrix<FloatType>::GetColumns() const {
  return num_columns_;
}

template <typename FloatType>
inline DynamicMatrix<FloatType>::DynamicMatrix(const DynamicMatrix<FloatType>& init) {
  this->num_columns_ = init.num_columns_;
  this->num_rows_ = init.num_rows_;
  theColumns = new DynamicVector<FloatType>[num_columns_];
  for (int i = 0; i < num_columns_; i++) {
    theColumns[i] = init.theColumns[i];
  }
}

template <typename FloatType>
inline DynamicVector<FloatType>& DynamicMatrix<FloatType>::operator[](const int& index) {
  return theColumns[index];
}

template <typename FloatType>
inline const DynamicVector<FloatType>& DynamicMatrix<FloatType>::operator[](const int& index) const {
  return theColumns[index];
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator+(const DynamicMatrix<FloatType>& op) const {
  DynamicMatrix<FloatType> result(num_columns_, num_rows_);
  for (int i = 0; i < num_columns_; i++) {
    result[i] = theColumns[i] + op.theColumns[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator-(const DynamicMatrix<FloatType>& op) const {
  DynamicMatrix<FloatType> result(num_columns_, num_rows_);
  for (int i = 0; i < num_columns_; i++) {
    result[i] = theColumns[i] - op.theColumns[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator-() const {
  DynamicMatrix<FloatType> result(num_columns_, num_rows_);
  for (int i = 0; i < num_columns_; i++) {
    result[i] = -theColumns[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator*(const FloatType& s) const {
  DynamicMatrix<FloatType> result(num_columns_, num_rows_);
  for (int i = 0; i < num_columns_; i++) {
    result[i] = theColumns[i] * s;
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator*(const DynamicMatrix<FloatType>& op) const {
  DynamicMatrix<FloatType> result(op.num_columns_, num_rows_);
  if (num_columns_ == op.num_rows_) {
    for (int r = 0; r < num_rows_; r++) {
      for (int c = 0; c < op.num_columns_; c++) {
        result[c][r] = 0.0;
        for (int i = 0; i < num_columns_; i++) {
          result[c][r] += (*this)[i][r] * op[c][i];
        }
      }
    }
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator/(const FloatType& s) const {
  DynamicMatrix<FloatType> result(num_columns_, num_rows_);
  for (int i = 0; i < num_columns_; i++) {
    result[i] = theColumns[i] / s;
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator+=(const DynamicMatrix<FloatType>& op) {
  for (int i = 0; i < num_columns_; i++) {
    theColumns[i] += op.theColumns[i];
  }
  return *this;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator-=(const DynamicMatrix<FloatType>& op) {
  for (int i = 0; i < num_columns_; i++) {
    theColumns[i] -= op.theColumns[i];
  }
  return *this;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator*=(const FloatType& op) {
  for (int i = 0; i < num_columns_; i++) {
    theColumns[i] *= op;
  }
  return *this;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator*=(const DynamicMatrix<FloatType>& op) {
  DynamicMatrix<FloatType> result = (*this) * op;
  *this = result;
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator/=(const FloatType& op) {

  for (int i = 0; i < num_columns_; i++) {
    theColumns[i] /= op;
  }
  return *this;
}

template <typename FloatType>
inline void DynamicMatrix<FloatType>::operator=(const DynamicMatrix<FloatType>& op) {
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    delete[] theColumns;
    theColumns = new DynamicVector<FloatType>[op.num_columns_];
  }
  this->num_columns_ = op.num_columns_;
  this->num_rows_ = op.num_rows_;
  for (int i = 0; i < op.num_columns_; i++) {
    theColumns[i] = op.theColumns[i];
  }
}

template <typename FloatType>
inline bool DynamicMatrix<FloatType>::operator==(const DynamicMatrix<FloatType>& op) const {
  for (int i = 0; i < num_columns_; i++) {
    if (theColumns[i] != op.theColumns[i]) return false;
  }
  return true;
}

template <typename FloatType>
inline bool DynamicMatrix<FloatType>::operator!=(const DynamicMatrix<FloatType>& op) const {
  for (int i = 0; i < num_columns_; i++) {
    if (theColumns[i] != op.theColumns[i]) return true;
  }
  return false;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicMatrix<FloatType>::operator*(const DynamicVector<FloatType>& v) const {
  DynamicVector<FloatType> result(num_rows_);
  for (int y = 0; y < num_rows_; y++) {
    result[y] = 0.0;
    for (int x = 0; x < num_columns_; x++) {
      result[y] += theColumns[x][y] * v[x];
    }
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::Transpose() const {
  DynamicMatrix<FloatType> result(num_rows_, num_columns_);
  for (int r = 0; r < num_rows_; r++) {
    for (int c = 0; c < num_columns_; c++) {
      result[r][c] = theColumns[c][r];
    }
  }
  return result;
}

template <typename FloatType>
inline FloatType DynamicMatrix<FloatType>::FrobeniusNorm() const {
  double sum_of_squred_norm = 0.0;
  for (int col_id = 0; col_id < num_columns_; col_id++) {
    DVectorD& buf = theColumns[col_id];
    sum_of_squred_norm += buf * buf;
  }
  return sqrt(sum_of_squred_norm);
}
// ---- destructor
template <typename FloatType>
inline DynamicMatrix<FloatType>::~DynamicMatrix() {
  delete[] theColumns;
}

typedef DynamicMatrix<float> DMatrixF;
typedef DynamicMatrix<double> DMatrixD;
typedef DynamicMatrix<int> DMatrixI;

#endif
