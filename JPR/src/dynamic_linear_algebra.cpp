
//---------------------------------------------------------------------------
#include "dynamic_linear_algebra.h"
//---------------------------------------------------------------------------

#include "math.h"
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

using namespace std;

template <typename FloatType>
inline DynamicVector<FloatType>::DynamicVector(int dim) {
  this->dim = dim;
  v = new FloatType[dim];
}

template <typename FloatType>
inline int DynamicVector<FloatType>::GetDim() const {
  return dim;
}

template <typename FloatType>
inline int DynamicVector<FloatType>::Size() const {
  return dim;
}

template <typename FloatType>
inline void DynamicVector<FloatType>::SetDim(const int dim) {
  this->dim = dim;
  delete[] v;
  v = new FloatType[dim];
}

template <typename FloatType>
inline void DynamicVector<FloatType>::SetZero() {
  for (int i = 0; i < dim; i++) {
    v[i] = (FloatType)0;
  }
}

template <typename FloatType>
inline DynamicVector<FloatType>::DynamicVector(const DynamicVector<FloatType>& init) {
  this->dim = init.dim;
  v = new FloatType[dim];
  for (int i = 0; i < dim; i++) {
    v[i] = init.v[i];
  }
}

template <typename FloatType>
inline FloatType& DynamicVector<FloatType>::operator[](const int& index) {
#ifdef DVECT_RANGE_CHECK
  if (index >= dim) {
    throw ERangeCheck("DynamicVector::operator[] - index exceeds bounds");
  }
#endif
  return v[index];
}

template <typename FloatType>
inline const FloatType& DynamicVector<FloatType>::operator[](const int& index) const {
#ifdef DVECT_RANGE_CHECK
  if (index >= dim) {
    throw ERangeCheck("DynamicVector::operator[] - index exceeds bounds");
  }
#endif
  return v[index];
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator+(const DynamicVector<FloatType>& op) const {
#ifdef DVECT_RANGE_CHECK
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator+ - incompatible dimensions.");
  }
#endif
  DynamicVector<FloatType> result(dim);
  for (int i = 0; i < dim; i++) {
    result[i] = v[i] + op.v[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator-(const DynamicVector<FloatType>& op) const {
#ifdef DVECT_RANGE_CHECK
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator- - incompatible dimensions.");
  }
#endif
  DynamicVector<FloatType> result(dim);
  for (int i = 0; i < dim; i++) {
    result[i] = v[i] - op.v[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator-() const {
  DynamicVector<FloatType> result(dim);
  for (int i = 0; i < dim; i++) {
    result[i] = -v[i];
  }
  return result;
}

template <typename FloatType>
inline FloatType DynamicVector<FloatType>::operator*(const DynamicVector<FloatType> op) const {
#ifdef DVECT_RANGE_CHECK
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator* - incompatible dimensions.");
  }
#endif
  FloatType result = 0.0;
  for (int i = 0; i < dim; i++) {
    result += v[i] * op.v[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator*(const FloatType& s) const {
  DynamicVector<FloatType> result(dim);
  for (int i = 0; i < dim; i++) {
    result[i] = s * v[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator/(const FloatType& s) const {
  DynamicVector<FloatType> result(dim);
  for (int i = 0; i < dim; i++) {
    result[i] = v[i] / s;
  }
  return result;
}

template <typename FloatType>
inline FloatType DynamicVector<FloatType>::ComponentSum() const {
  FloatType result = 0;
  for (int i = 0; i < dim; i++) {
    result += v[i];
  }
  return result;
}

template <typename FloatType>
inline FloatType DynamicVector<FloatType>::Average() const {
  return ComponentSum() / dim;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator+=(const DynamicVector<FloatType>& op) {
#ifdef DVECT_RANGE_CHECK
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator+= - incompatible dimensions.");
  }
#endif
  for (int i = 0; i < dim; i++) {
    v[i] += op.v[i];
  }
  return *this;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator-=(const DynamicVector<FloatType>& op) {
#ifdef DVECT_RANGE_CHECK
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator-= - incompatible dimensions.");
  }
#endif
  for (int i = 0; i < dim; i++) {
    v[i] -= op.v[i];
  }
  return *this;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator*=(const FloatType& s) {
  for (int i = 0; i < dim; i++) {
    v[i] *= s;
  }
  return *this;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::operator/=(const FloatType& s) {
  for (int i = 0; i < dim; i++) {
    v[i] /= s;
  }
  return *this;
}

template <typename FloatType>
inline void DynamicVector<FloatType>::operator=(const DynamicVector<FloatType>& op) {
  if (dim != op.dim) {
    delete[] v;
    v = new FloatType[op.dim];
  }
  this->dim = op.dim;
  for (int i = 0; i < op.dim; i++) {
    v[i] = op.v[i];
  }
}

template <typename FloatType>
inline bool DynamicVector<FloatType>::operator==(const DynamicVector<FloatType> op) const {
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator== - incompatible dimensions.");
  }
  for (int i = 0; i < dim; i++) {
    if (v[i] != op.v[i]) return false;
  }
  return true;
}

template <typename FloatType>
inline bool DynamicVector<FloatType>::operator!=(const DynamicVector<FloatType> op) const {
#ifdef DVECT_RANGE_CHECK
  if (dim != op.dim) {
    throw ERangeCheck("DynamicVector::operator!= - incompatible dimensions.");
  }
#endif
  for (int i = 0; i < dim; i++) {
    if (v[i] != op.v[i]) return true;
  }
  return false;
}

template <typename FloatType>
inline FloatType DynamicVector<FloatType>::MaxComponent() const {
  FloatType result = v[0];
  for (int i = 1; i < dim; i++) {
    if (v[i] > result) {
      result = v[i];
    }
  }
  return result;
}

template <typename FloatType>
inline FloatType DynamicVector<FloatType>::MinComponent() const {
  FloatType result = v[0];
  for (int i = 1; i < dim; i++) {
    if (v[i] < result) {
      result = v[i];
    }
  }
  return result;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicVector<FloatType>::AbsoluteValues() const {
  DynamicVector<FloatType> result(dim);
  for (int i = 0; i < dim; i++) {
    result[i] = fabs(v[i]);
  }
  return result;
}

template <typename FloatType>
inline DynamicVector<FloatType>::~DynamicVector() {
  delete[] v;
}

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
#ifdef DVECT_RANGE_CHECK
  if (index >= num_columns_) {
    throw ERangeCheck("DynamicMatrix::operator[] - index exceeds bounds");
  }
#endif
  return theColumns[index];
}

template <typename FloatType>
inline const DynamicVector<FloatType>& DynamicMatrix<FloatType>::operator[](const int& index) const {
#ifdef DVECT_RANGE_CHECK
  if (index >= num_columns_) {
    throw ERangeCheck("DynamicMatrix::operator[] - index exceeds bounds");
  }
#endif
  return theColumns[index];
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator+(const DynamicMatrix<FloatType>& op) const {
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator+ - incompatible dimensions.");
  }
#endif
  DynamicMatrix<FloatType> result(num_columns_, num_rows_);
  for (int i = 0; i < num_columns_; i++) {
    result[i] = theColumns[i] + op.theColumns[i];
  }
  return result;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator-(const DynamicMatrix<FloatType>& op) const {
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator- - incompatible dimensions.");
  }
#endif
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
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator+= - incompatible dimensions.");
  }
#endif
  for (int i = 0; i < num_columns_; i++) {
    theColumns[i] += op.theColumns[i];
  }
  return *this;
}

template <typename FloatType>
inline DynamicMatrix<FloatType> DynamicMatrix<FloatType>::operator-=(const DynamicMatrix<FloatType>& op) {
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator-= - incompatible dimensions.");
  }
#endif
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
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator*= - incompatible dimensions.");
  }
#endif
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
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator== - incompatible dimensions.");
  }
#endif
  for (int i = 0; i < num_columns_; i++) {
    if (theColumns[i] != op.theColumns[i]) return false;
  }
  return true;
}

template <typename FloatType>
inline bool DynamicMatrix<FloatType>::operator!=(const DynamicMatrix<FloatType>& op) const {
#ifdef DVECT_RANGE_CHECK
  if (num_columns_ != op.num_columns_ || num_rows_ != op.num_rows_) {
    throw ERangeCheck("DynamicMatrix::operator!= - incompatible dimensions.");
  }
#endif
  for (int i = 0; i < num_columns_; i++) {
    if (theColumns[i] != op.theColumns[i]) return true;
  }
  return false;
}

template <typename FloatType>
inline DynamicVector<FloatType> DynamicMatrix<FloatType>::operator*(const DynamicVector<FloatType>& v) const {
#ifdef DVECT_RANGE_CHECK
  if (v.getDim() != num_columns_) {
    throw ERangeCheck("DynamicMatrix::operator* (vect) - incompatible dimensions.");
  }
#endif
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