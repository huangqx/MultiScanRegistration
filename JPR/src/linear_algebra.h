//---------------------------------------------------------------------------
#ifndef linear_algebra_h_
#define linear_algebra_h_

#include "math.h"

template <class FloatType, unsigned dim>
class StaticVector {
 private:
   FloatType v[dim];
 public:
   inline FloatType& operator[](const unsigned &index);
   inline const FloatType& operator[](const unsigned &index) const;

   //  ----  operators +, -, *, /, crossProduct
   inline StaticVector<FloatType, dim> operator+(const StaticVector<FloatType, dim> &op) const;
   inline StaticVector<FloatType, dim> operator-(const StaticVector<FloatType, dim> &op) const;
   inline StaticVector<FloatType, dim> operator-() const;
   inline FloatType operator*(const StaticVector<FloatType, dim> op) const;
   inline StaticVector<FloatType, dim> operator*(const FloatType &s) const;
   inline StaticVector<FloatType, dim> operator/(const FloatType &s) const;
   StaticVector<FloatType, dim> crossProduct(const StaticVector<FloatType, dim> &op) const;
   inline StaticVector<FloatType, dim> componentProduct(const StaticVector<FloatType, dim> &op) const;

   //  ---- operators +=, -=, *=, /=
   inline StaticVector<FloatType, dim> operator+=(const StaticVector<FloatType, dim> &op);
   inline StaticVector<FloatType, dim> operator-=(const StaticVector<FloatType, dim> &op);
   inline StaticVector<FloatType, dim> operator*=(const FloatType &s);
   inline StaticVector<FloatType, dim> operator/=(const FloatType &s);

   //  ---- operators =, ==, !=
   inline bool operator==(const StaticVector<FloatType, dim> op) const;
   inline bool operator!=(const StaticVector<FloatType, dim> op) const;

   inline FloatType* data();
   inline const FloatType* data() const;
   inline void normalize();
   inline FloatType getSqrNorm() const;

   /// compatibility with DynamicVector
   static inline unsigned getDim() {return dim;}
   static inline unsigned size() {return dim;}
};


// ---------------------------------------------------------------------------
//
//                                 StaticMatrix
//
// ---------------------------------------------------------------------------


template <class FloatType, unsigned columns, unsigned rows>
class StaticMatrix {
 private:
   StaticVector<FloatType, rows> theColumns[columns];

 public:
   // ---- constructors
   inline StaticMatrix<FloatType, columns, rows>();

   // ---- element access
   inline StaticVector<FloatType, rows>& operator[](const unsigned &index);
   inline const StaticVector<FloatType, rows>& operator[](const unsigned &index) const;

   //  ----  operators +, -, *, /
   inline StaticMatrix<FloatType, columns, rows> operator+(const StaticMatrix<FloatType, columns, rows> &op) const;
   inline StaticMatrix<FloatType, columns, rows> operator-(const StaticMatrix<FloatType, columns, rows> &op) const;
   inline StaticMatrix<FloatType, columns, rows> operator-() const;
   inline StaticMatrix<FloatType, columns, rows> operator*(const FloatType &s) const;
   /* This:    colums: columns (b)
               rows: rows (a)
      Operand: colums: colums2 (c)
               rows: columns (b)
      Result:  colums: columns2 (a)
               rows: rows (c)
   */
   template <unsigned columns2>
   inline StaticMatrix<FloatType, columns2, rows> operator*( const StaticMatrix<FloatType, columns2, columns>& mat ) const;
   
   inline StaticMatrix<FloatType, columns, rows> operator/(const FloatType &s) const;

   //  ---- operators +=, -=, *=, /=
   inline StaticMatrix<FloatType, columns, rows> operator+=(const StaticMatrix<FloatType, columns, rows> &op);
   inline StaticMatrix<FloatType, columns, rows> operator-=(const StaticMatrix<FloatType, columns, rows> &op);
   inline StaticMatrix<FloatType, columns, rows> operator*=(const FloatType &op);
   inline StaticMatrix<FloatType, columns, rows> operator*=(const StaticMatrix<FloatType, columns, rows> &op);
   inline StaticMatrix<FloatType, columns, rows> operator/=(const FloatType &op);

   //  ---- operators =, ==, !=
   inline StaticMatrix<FloatType, columns, rows> operator=(const StaticMatrix<FloatType, columns, rows> &op);
   inline bool operator==(const StaticMatrix<FloatType, columns, rows> &op) const;
   inline bool operator!=(const StaticMatrix<FloatType, columns, rows> &op) const;

   // ---- multiplying with vectors
   inline StaticVector<FloatType, rows> operator*(const StaticVector<FloatType, columns> &v) const;

   inline FloatType* data();
   inline const FloatType* data() const;
   inline void changeRows(unsigned row1, unsigned row2);
   inline void multRow(unsigned row, FloatType value);
   inline void combineRows(unsigned row, unsigned with, FloatType by);
   inline void addRows(unsigned row, unsigned with);
   inline StaticMatrix<FloatType, rows, columns> transpose() const;
   inline FloatType getDeterminant() const;

   /// Sort eigenvectors according to eigenvalues in descending order
   void eigenSort( StaticVector<FloatType,rows>& eigenValues, StaticMatrix<FloatType,rows,columns>& eigenVectors ) const;
   
   /// compatibility with DynamicVector
   static inline unsigned getRowsDim() {return rows;}
   static inline unsigned getColsDim() {return columns;}
   static inline unsigned getRows()    {return rows;}
   static inline unsigned getColumns() {return columns;}
};


// ---------------------------------------------------------------------------
//
//                                  instances
//
// ---------------------------------------------------------------------------


typedef StaticVector<float,2> Vector2f;
typedef StaticVector<float,3> Vector3f;
typedef StaticVector<float,4> Vector4f;
typedef StaticVector<float,5> Vector5f;
typedef StaticVector<float,6> Vector6f;
typedef StaticVector<float,20> Vector20f;
typedef StaticVector<double,2> Vector2d;
typedef StaticVector<double,3> Vector3d;
typedef StaticVector<double,4> Vector4d;
typedef StaticVector<double,5> Vector5d;
typedef StaticVector<double,6> Vector6d;
typedef StaticVector<double,8> Vector8d;
typedef StaticVector<double,12> Vector12d;
typedef StaticVector<double,20> Vector20d;

typedef StaticMatrix<float,2,2> Matrix2f;
typedef StaticMatrix<float,3,3> Matrix3f;
typedef StaticMatrix<float,4,4> Matrix4f;
typedef StaticMatrix<float,5,5> Matrix5f;
typedef StaticMatrix<double,2,2> Matrix2d;
typedef StaticMatrix<double,3,3> Matrix3d;
typedef StaticMatrix<double,4,4> Matrix4d;
typedef StaticMatrix<double,5,5> Matrix5d;
typedef StaticMatrix<double,6,6> Matrix6d;
typedef StaticMatrix<double,8,8> Matrix8d;
typedef StaticMatrix<double,12,12> Matrix12d;
typedef StaticMatrix<double,6,12> Matrix12x6d;
typedef StaticMatrix<double,12,6> Matrix6x12d;


// ---------------------------------------------------------------------------
//
//                                 misc. functions
//
// ---------------------------------------------------------------------------

inline void CompleteCoordinateFrame(const Vector3f &axis_z,
                             Vector3f *axis_x,
                             Vector3f *axis_y);

inline float TriangleArea(const Vector3f &p1,
                          const Vector3f &p2,
                          const Vector3f &p3);

inline double VectorAngle(const Vector3f &vec, 
                   const Vector3f &axis_x,
                   const Vector3f &axis_y);

inline double VectorAngle(const Vector3d &vec, 
                          const Vector3d &axis_x,
                          const Vector3d &axis_y);

inline double AngleBetweenVectors(const Vector3f &origin,
                                  const Vector3f &v1,
                                  const Vector3f &v2);

inline double AngleBetweenVectors(const Vector3d &origin,
                                  const Vector3d &v1,
                                  const Vector3d &v2);


template <class FloatType, unsigned dim>
inline FloatType& StaticVector<FloatType, dim>::operator[](const unsigned& index) {
#ifdef RANGE_CHECK
  if (index >= dim) {
    throw ERangeCheck("StaticVector::operator[] - index exceeds bounds");
  }
#endif
  return v[index];
}

template <class FloatType, unsigned dim>
inline const FloatType& StaticVector<FloatType, dim>::operator[](const unsigned& index) const {
#ifdef RANGE_CHECK
  if (index >= dim) {
    throw ERangeCheck("StaticVector::operator[] - index exceeds bounds");
  }
#endif
  return v[index];
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator+(const StaticVector<FloatType, dim>& op) const {
  StaticVector<FloatType, dim> result;
  for (unsigned i = 0; i < dim; i++) {
    result[i] = v[i] + op.v[i];
  }
  return result;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator-(const StaticVector<FloatType, dim>& op) const {
  StaticVector<FloatType, dim> result;
  for (unsigned i = 0; i < dim; i++) {
    result[i] = v[i] - op.v[i];
  }
  return result;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator-() const {
  StaticVector<FloatType, dim> result;
  for (unsigned i = 0; i < dim; i++) {
    result[i] = -v[i];
  }
  return result;
}

template <class FloatType, unsigned dim>
inline FloatType StaticVector<FloatType, dim>::operator*(const StaticVector<FloatType, dim> op) const {
  FloatType result = 0.0;
  for (unsigned i = 0; i < dim; i++) {
    result += v[i] * op.v[i];
  }
  return result;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator*(const FloatType& s) const {
  StaticVector<FloatType, dim> result;
  for (unsigned i = 0; i < dim; i++) {
    result[i] = s * v[i];
  }
  return result;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator/(const FloatType& s) const {
  StaticVector<FloatType, dim> result;
  for (unsigned i = 0; i < dim; i++) {
    result[i] = v[i] / s;
  }
  return result;
}

template <class FloatType, unsigned dim>
StaticVector<FloatType, dim> StaticVector<FloatType, dim>::crossProduct(const StaticVector<FloatType, dim>& op) const {
  StaticVector<FloatType, dim> result;
  result[0] = v[1] * op.v[2] - v[2] * op.v[1];
  result[1] = v[2] * op.v[0] - v[0] * op.v[2];
  result[2] = v[0] * op.v[1] - v[1] * op.v[0];
  return result;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::componentProduct(const StaticVector<FloatType, dim>& op) const {
  StaticVector<FloatType, dim> result;
  for (unsigned i = 0; i < dim; i++) {
    result[i] = v[i] * op.v[i];
  }
  return result;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator+=(const StaticVector<FloatType, dim>& op) {
  for (unsigned i = 0; i < dim; i++) {
    v[i] += op.v[i];
  }
  return *this;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator-=(const StaticVector<FloatType, dim>& op) {
  for (unsigned i = 0; i < dim; i++) {
    v[i] -= op.v[i];
  }
  return *this;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator*=(const FloatType& s) {
  for (unsigned i = 0; i < dim; i++) {
    v[i] *= s;
  }
  return *this;
}

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> StaticVector<FloatType, dim>::operator/=(const FloatType& s) {
  for (unsigned i = 0; i < dim; i++) {
    v[i] /= s;
  }
  return *this;
}

template <class FloatType, unsigned dim>
inline bool StaticVector<FloatType, dim>::operator==(const StaticVector<FloatType, dim> op) const {
  for (unsigned i = 0; i < dim; i++) {
    if (v[i] != op.v[i]) return false;
  }
  return true;
}

template <class FloatType, unsigned dim>
inline bool StaticVector<FloatType, dim>::operator!=(const StaticVector<FloatType, dim> op) const {
  for (unsigned i = 0; i < dim; i++) {
    if (v[i] != op.v[i]) return true;
  }
  return false;
}

template <class FloatType, unsigned dim>
inline FloatType* StaticVector<FloatType, dim>::data() {
  return &(v[0]);
}

template <class FloatType, unsigned dim>
inline const FloatType* StaticVector<FloatType, dim>::data() const {
  return &(v[0]);
}

template <class FloatType, unsigned dim>
inline void StaticVector<FloatType, dim>::normalize() {
  FloatType norm = 0.0;
  for (int i = 0; i < dim; i++) {
    norm += v[i] * v[i];
  }
  if (norm > 1e-12f)
    (*this) = (*this) / static_cast<float> (sqrt(norm));
}

template <class FloatType, unsigned dim>
inline FloatType StaticVector<FloatType, dim>::getSqrNorm() const {
  FloatType norm = (FloatType)0.0;
  for (int i = 0; i < dim; i++) {
    norm += v[i] * v[i];
  }
  return norm;
}

template<>
inline StaticVector<float, 3> StaticVector<float, 3>::crossProduct(const StaticVector<float, 3>& op) const {
  StaticVector<float, 3> result;
  result.v[0] = v[1] * op.v[2] - v[2] * op.v[1];
  result.v[1] = v[2] * op.v[0] - v[0] * op.v[2];
  result.v[2] = v[0] * op.v[1] - v[1] * op.v[0];
  return result;
}

template<>
inline StaticVector<int, 3> StaticVector<int, 3>::crossProduct(const StaticVector<int, 3>& op) const {
  StaticVector<int, 3> result;
  result.v[0] = v[1] * op.v[2] - v[2] * op.v[1];
  result.v[1] = v[2] * op.v[0] - v[0] * op.v[2];
  result.v[2] = v[0] * op.v[1] - v[1] * op.v[0];
  return result;
}


// ---- helper functions


template <class FloatTypeDest, class FloatTypeSource, int Dim>
StaticVector<FloatTypeDest, Dim> convertVector(const StaticVector<FloatTypeSource, Dim>& source) {
  StaticVector<FloatTypeDest, Dim> result;
  for (int i = 0; i < Dim; i++) {
    result[i] = (FloatTypeDest)source[i];
  }
  return result;
}



// ---------------------------------------------------------------------------
//
//                                 StaticMatrix
//
// ---------------------------------------------------------------------------



template <class FloatType, unsigned columns, unsigned rows>
inline StaticVector<FloatType, rows>& StaticMatrix<FloatType, columns, rows>::operator[](const unsigned& index) {
#ifdef RANGE_CHECK
  if (index >= columns) {
    throw ERangeCheck("StaticMatrix::operator[] - index exceeds bounds");
  }
#endif
  return theColumns[index];
}

template <class FloatType, unsigned columns, unsigned rows>
inline const StaticVector<FloatType, rows>& StaticMatrix<FloatType, columns, rows>::operator[](const unsigned& index) const {
#ifdef RANGE_CHECK
  if (index >= columns) {
    throw ERangeCheck("StaticMatrix::operator[] - index exceeds bounds");
  }
#endif
  return theColumns[index];
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator+(const StaticMatrix<FloatType, columns, rows>& op) const {
  StaticMatrix<FloatType, columns, rows> result;
  for (unsigned i = 0; i < columns; i++) {
    result[i] = theColumns[i] + op.theColumns[i];
  }
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator-(const StaticMatrix<FloatType, columns, rows>& op) const {
  StaticMatrix<FloatType, columns, rows> result;
  for (unsigned i = 0; i < columns; i++) {
    result[i] = theColumns[i] - op.theColumns[i];
  }
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator-() const {
  StaticMatrix<FloatType, columns, rows> result;
  for (unsigned i = 0; i < columns; i++) {
    result[i] = -theColumns[i];
  }
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator*(const FloatType& s) const {
  StaticMatrix<FloatType, columns, rows> result;
  for (unsigned i = 0; i < columns; i++) {
    result[i] = theColumns[i] * s;
  }
  return result;
}


/* This:    colums: columns
            rows: rows
   Operand: colums: colums2
            rows: columns
   Result:  colums: columns2
            rows: rows
*/
template <class FloatType, unsigned columns, unsigned rows>
template <unsigned columns2>
StaticMatrix<FloatType, columns2, rows> StaticMatrix<FloatType, columns, rows>::operator*(const StaticMatrix<FloatType, columns2, columns>& mat) const
{
  StaticMatrix<FloatType, columns2, rows> res;
  for (unsigned char row = 0; row < rows; row++) {
    for (unsigned char column = 0; column < columns2; column++) {
      res[column][row] = 0;
      for (unsigned char k = 0; k < columns; k++)
        res[column][row] += theColumns[k][row] * mat[column][k];
    }
  }
  return res;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator/(const FloatType& s) const {
  StaticMatrix<FloatType, columns, rows> result;
  for (unsigned i = 0; i < columns; i++) {
    result[i] = theColumns[i] / s;
  }
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator+=(const StaticMatrix<FloatType, columns, rows>& op) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i] += op.theColumns[i];
  }
  return *this;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator-=(const StaticMatrix<FloatType, columns, rows>& op) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i] -= op.theColumns[i];
  }
  return *this;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator*=(const FloatType& op) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i] *= op;
  }
  return *this;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator*=(const StaticMatrix<FloatType, columns, rows>& op) {
  StaticMatrix<FloatType, columns, rows> result = (*this) * op;
  *this = result;
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator/=(const FloatType& op) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i] /= op;
  }
  return *this;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows> StaticMatrix<FloatType, columns, rows>::operator=(const StaticMatrix<FloatType, columns, rows>& op) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i] = op.theColumns[i];
  }
  return *this;
}

template <class FloatType, unsigned columns, unsigned rows>
inline bool StaticMatrix<FloatType, columns, rows>::operator==(const StaticMatrix<FloatType, columns, rows>& op) const {
  for (unsigned i = 0; i < columns; i++) {
    if (theColumns[i] != op.theColumns[i]) return false;
  }
  return true;
}

template <class FloatType, unsigned columns, unsigned rows>
inline bool StaticMatrix<FloatType, columns, rows>::operator!=(const StaticMatrix<FloatType, columns, rows>& op) const {
  for (unsigned i = 0; i < columns; i++) {
    if (theColumns[i] != op.theColumns[i]) return true;
  }
  return false;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticVector<FloatType, rows> StaticMatrix<FloatType, columns, rows>::operator*(const StaticVector<FloatType, columns>& v) const {
  StaticVector<FloatType, rows> result;
  for (unsigned y = 0; y < rows; y++) {
    result[y] = 0.0;
    for (unsigned x = 0; x < columns; x++) {
      result[y] += theColumns[x][y] * v[x];
    }
  }
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, columns, rows>::StaticMatrix() {
  for (unsigned j = 0; j < columns; j++) {
    for (unsigned i = 0; i < rows; i++) {
      if (i == j)
        theColumns[j][i] = 1.0;
      else
        theColumns[j][i] = 0.0;
    }
  }
}


template <class FloatType, unsigned columns, unsigned rows>
inline FloatType* StaticMatrix<FloatType, columns, rows>::data() {
  return &(theColumns[0][0]);
}

template <class FloatType, unsigned columns, unsigned rows>
inline const FloatType* StaticMatrix<FloatType, columns, rows>::data() const {
  return &(theColumns[0][0]);
}

template <class FloatType, unsigned columns, unsigned rows>
inline void StaticMatrix<FloatType, columns, rows>::changeRows(unsigned row1, unsigned row2) {
  StaticVector<FloatType, columns> h;
  unsigned i;
  for (i = 0; i < columns; i++) {
    h[i] = theColumns[i][row2];
  }
  for (i = 0; i < columns; i++) {
    theColumns[i][row2] = theColumns[i][row1];
  }
  for (i = 0; i < columns; i++) {
    theColumns[i][row1] = h[i];
  }
}

template <class FloatType, unsigned columns, unsigned rows>
inline void StaticMatrix<FloatType, columns, rows>::multRow(unsigned row, FloatType value) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i][row] *= value;
  }
}

template <class FloatType, unsigned columns, unsigned rows>
inline void StaticMatrix<FloatType, columns, rows>::combineRows(unsigned row, unsigned with, FloatType by) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i][row] += theColumns[i][with] * by;
  }
}

template <class FloatType, unsigned columns, unsigned rows>
inline void StaticMatrix<FloatType, columns, rows>::addRows(unsigned row, unsigned with) {
  for (unsigned i = 0; i < columns; i++) {
    theColumns[i][row] += theColumns[i][with];
  }
}

template <class FloatType, unsigned columns, unsigned rows>
inline StaticMatrix<FloatType, rows, columns> StaticMatrix<FloatType, columns, rows>::transpose() const {
  StaticMatrix<FloatType, rows, columns> result;
  for (unsigned r = 0; r < rows; r++) {
    for (unsigned c = 0; c < columns; c++) {
      result[r][c] = theColumns[c][r];
    }
  }
  return result;
}

template <class FloatType, unsigned columns, unsigned rows>
inline FloatType StaticMatrix<FloatType, columns, rows>::getDeterminant() const {
  if ((rows == 2) && (columns == 2))
    return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
  else if ((rows == 3) && (columns == 3))
    return (*this)[0][0] * ((*this)[1][1] * (*this)[2][2] - (*this)[1][2] * (*this)[2][1]) -
    (*this)[0][1] * ((*this)[1][0] * (*this)[2][2] - (*this)[1][2] * (*this)[2][0]) +
    (*this)[0][2] * ((*this)[1][0] * (*this)[2][1] - (*this)[1][1] * (*this)[2][0]);
  else if (rows == 1 && columns == 1)
    return (*this)[0][0];
  else
    // limitiation: does not work for dimensions higher larger 3
  {

  }
}

/**
 *  Helping structure for jacobi iteration
 */
#define STATIC_MATRIX_ROTATE(a,i,j,k,l) g=a[i][j];\
                                        h=a[k][l];\
                                        a[i][j]=g-s*(h+g*tau);\
                                        a[k][l]=h+s*(g-h*tau);

template <class FloatType, unsigned columns, unsigned rows>
void StaticMatrix<FloatType, columns, rows>::eigenSort(StaticVector<FloatType, rows>& eigenValues, StaticMatrix<FloatType, rows, columns>& eigenVectors) const
{
  int i, j, k;
  FloatType p;
  for (i = 0; i < rows - 1; i++) {
    p = eigenValues[k = i];
    for (j = i + 1; j < rows; j++)
      if (eigenValues[j] >= p) p = eigenValues[k = j];
    if (k != i) {
      eigenValues[k] = eigenValues[i];
      eigenValues[i] = p;
      for (j = 0; j < rows; j++) {
        p = eigenVectors[j][i];
        eigenVectors[j][i] = eigenVectors[j][k];
        eigenVectors[j][k] = p;
      }
    }
  }
}


// --- explicit specialization for dim 2,3,4 cases:


#ifdef LINEAR_ALGEBRA_EXPLICIT_LOOP_UNROOLING


template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator+(const StaticVector<float, 4>& op) const {
  StaticVector<float, 4> result;
  result.v[0] = v[0] + op.v[0];
  result.v[1] = v[1] + op.v[1];
  result.v[2] = v[2] + op.v[2];
  result.v[3] = v[3] + op.v[3];
  return result;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator+(const StaticVector<float, 3>& op) const {
  StaticVector<float, 3> result;
  result.v[0] = v[0] + op.v[0];
  result.v[1] = v[1] + op.v[1];
  result.v[2] = v[2] + op.v[2];
  return result;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator+(const StaticVector<float, 2>& op) const {
  StaticVector<float, 2> result;
  result.v[0] = v[0] + op.v[0];
  result.v[1] = v[1] + op.v[1];
  return result;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator-(const StaticVector<float, 4>& op) const {
  StaticVector<float, 4> result;
  result.v[0] = v[0] - op.v[0];
  result.v[1] = v[1] - op.v[1];
  result.v[2] = v[2] - op.v[2];
  result.v[3] = v[3] - op.v[3];
  return result;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator-(const StaticVector<float, 3>& op) const {
  StaticVector<float, 3> result;
  result.v[0] = v[0] - op.v[0];
  result.v[1] = v[1] - op.v[1];
  result.v[2] = v[2] - op.v[2];
  return result;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator-(const StaticVector<float, 2>& op) const {
  StaticVector<float, 2> result;
  result.v[0] = v[0] - op.v[0];
  result.v[1] = v[1] - op.v[1];
  return result;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator-() const {
  StaticVector<float, 4> result;
  result.v[0] = -v[0];
  result.v[1] = -v[1];
  result.v[2] = -v[2];
  result.v[3] = -v[3];
  return result;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator-() const {
  StaticVector<float, 3> result;
  result.v[0] = -v[0];
  result.v[1] = -v[1];
  result.v[2] = -v[2];
  return result;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator-() const {
  StaticVector<float, 2> result;
  result.v[0] = -v[0];
  result.v[1] = -v[1];
  return result;
}

template<>
inline float StaticVector<float, 4>::operator*(const StaticVector<float, 4> op) const {
  float result = 0.0;
  result += v[0] * op.v[0];
  result += v[1] * op.v[1];
  result += v[2] * op.v[2];
  result += v[3] * op.v[3];
  return result;
}
template<>
inline float StaticVector<float, 3>::operator*(const StaticVector<float, 3> op) const {
  float result = 0.0;
  result += v[0] * op.v[0];
  result += v[1] * op.v[1];
  result += v[2] * op.v[2];
  return result;
}
template<>
inline float StaticVector<float, 2>::operator*(const StaticVector<float, 2> op) const {
  float result = 0.0;
  result += v[0] * op.v[0];
  result += v[1] * op.v[1];
  return result;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator*(const float& s) const {
  StaticVector<float, 4> result;
  result.v[0] = s * v[0];
  result.v[1] = s * v[1];
  result.v[2] = s * v[2];
  result.v[3] = s * v[3];
  return result;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator*(const float& s) const {
  StaticVector<float, 3> result;
  result.v[0] = s * v[0];
  result.v[1] = s * v[1];
  result.v[2] = s * v[2];
  return result;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator*(const float& s) const {
  StaticVector<float, 2> result;
  result.v[0] = s * v[0];
  result.v[1] = s * v[1];
  return result;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator/(const float& s) const {
  StaticVector<float, 4> result;
  result.v[0] = v[0] / s;
  result.v[1] = v[1] / s;
  result.v[2] = v[2] / s;
  result.v[3] = v[3] / s;
  return result;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator/(const float& s) const {
  StaticVector<float, 3> result;
  result.v[0] = v[0] / s;
  result.v[1] = v[1] / s;
  result.v[2] = v[2] / s;
  return result;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator/(const float& s) const {
  StaticVector<float, 2> result;
  result.v[0] = v[0] / s;
  result.v[1] = v[1] / s;
  return result;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::componentProduct(const StaticVector<float, 4>& op) const {
  StaticVector<float, 4> result;
  result.v[0] = v[0] * op.v[0];
  result.v[1] = v[1] * op.v[1];
  result.v[2] = v[2] * op.v[2];
  result.v[3] = v[3] * op.v[3];
  return result;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::componentProduct(const StaticVector<float, 3>& op) const {
  StaticVector<float, 3> result;
  result.v[0] = v[0] * op.v[0];
  result.v[1] = v[1] * op.v[1];
  result.v[2] = v[2] * op.v[2];
  return result;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::componentProduct(const StaticVector<float, 2>& op) const {
  StaticVector<float, 2> result;
  result.v[0] = v[0] * op.v[0];
  result.v[1] = v[1] * op.v[1];
  return result;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator+=(const StaticVector<float, 4>& op) {
  v[0] += op.v[0];
  v[1] += op.v[1];
  v[2] += op.v[2];
  v[3] += op.v[3];
  return *this;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator+=(const StaticVector<float, 3>& op) {
  v[0] += op.v[0];
  v[1] += op.v[1];
  v[2] += op.v[2];
  return *this;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator+=(const StaticVector<float, 2>& op) {
  v[0] += op.v[0];
  v[1] += op.v[1];
  return *this;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator-=(const StaticVector<float, 4>& op) {
  v[0] -= op.v[0];
  v[1] -= op.v[1];
  v[2] -= op.v[2];
  v[3] -= op.v[3];
  return *this;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator-=(const StaticVector<float, 3>& op) {
  v[0] -= op.v[0];
  v[1] -= op.v[1];
  v[2] -= op.v[2];
  return *this;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator-=(const StaticVector<float, 2>& op) {
  v[0] -= op.v[0];
  v[1] -= op.v[1];
  return *this;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator*=(const float& s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
  v[3] *= s;
  return *this;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator*=(const float& s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
  return *this;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator*=(const float& s) {
  v[0] *= s;
  v[1] *= s;
  return *this;
}

template<>
inline StaticVector<float, 4> StaticVector<float, 4>::operator/=(const float& s) {
  v[0] /= s;
  v[1] /= s;
  v[2] /= s;
  v[3] /= s;
  return *this;
}
template<>
inline StaticVector<float, 3> StaticVector<float, 3>::operator/=(const float& s) {
  v[0] /= s;
  v[1] /= s;
  v[2] /= s;
  return *this;
}
template<>
inline StaticVector<float, 2> StaticVector<float, 2>::operator/=(const float& s) {
  v[0] /= s;
  v[1] /= s;
  return *this;
}

template<>
inline bool StaticVector<float, 4>::operator==(const StaticVector<float, 4> op) const {
  if (v[0] != op.v[0]) return false;
  if (v[1] != op.v[1]) return false;
  if (v[2] != op.v[2]) return false;
  if (v[3] != op.v[3]) return false;
  return true;
}
template<>
inline bool StaticVector<float, 3>::operator==(const StaticVector<float, 3> op) const {
  if (v[0] != op.v[0]) return false;
  if (v[1] != op.v[1]) return false;
  if (v[2] != op.v[2]) return false;
  return true;
}
template<>
inline bool StaticVector<float, 2>::operator==(const StaticVector<float, 2> op) const {
  if (v[0] != op.v[0]) return false;
  if (v[1] != op.v[1]) return false;
  return true;
}

template<>
inline bool StaticVector<float, 4>::operator!=(const StaticVector<float, 4> op) const {
  if (v[0] != op.v[0]) return true;
  if (v[1] != op.v[1]) return true;
  if (v[2] != op.v[2]) return true;
  if (v[3] != op.v[3]) return true;
  return false;
}
template<>
inline bool StaticVector<float, 3>::operator!=(const StaticVector<float, 3> op) const {
  if (v[0] != op.v[0]) return true;
  if (v[1] != op.v[1]) return true;
  if (v[2] != op.v[2]) return true;
  return false;
}
template<>
inline bool StaticVector<float, 2>::operator!=(const StaticVector<float, 2> op) const {
  if (v[0] != op.v[0]) return true;
  if (v[1] != op.v[1]) return true;
  return false;
}



template<>
inline StaticMatrix<float, 2, 2>::StaticMatrix() {
  theColumns[0][0] = 1.0f; theColumns[0][1] = 0.0f;
  theColumns[1][0] = 0.0f; theColumns[1][1] = 1.0f;
}

template<>
inline StaticMatrix<float, 3, 3>::StaticMatrix() {
  theColumns[0][0] = 1.0f; theColumns[0][1] = 0.0f; theColumns[0][2] = 0.0f;
  theColumns[1][0] = 0.0f; theColumns[1][1] = 1.0f; theColumns[1][2] = 0.0f;
  theColumns[2][0] = 0.0f; theColumns[2][1] = 0.0f; theColumns[2][2] = 1.0f;
}

template<>
inline StaticMatrix<float, 4, 4>::StaticMatrix() {
  theColumns[0][0] = 1.0f; theColumns[0][1] = 0.0f; theColumns[0][2] = 0.0f; theColumns[0][3] = 0.0f;
  theColumns[1][0] = 0.0f; theColumns[1][1] = 1.0f; theColumns[1][2] = 0.0f; theColumns[1][3] = 0.0f;
  theColumns[2][0] = 0.0f; theColumns[2][1] = 0.0f; theColumns[2][2] = 1.0f; theColumns[2][3] = 0.0f;
  theColumns[3][0] = 0.0f; theColumns[3][1] = 0.0f; theColumns[3][2] = 0.0f; theColumns[3][3] = 1.0f;
}

template<>
inline StaticMatrix<float, 2, 3>::StaticMatrix() {
  theColumns[0][0] = 0.0f; theColumns[0][1] = 0.0f;
  theColumns[1][0] = 0.0f; theColumns[1][1] = 0.0f;
  theColumns[2][0] = 0.0f; theColumns[2][1] = 0.0f;
}

template<>
inline StaticMatrix<float, 2, 2> StaticMatrix<float, 2, 2>::operator+(const StaticMatrix<float, 2, 2>& op) const {
  StaticMatrix<float, 2, 2> result;
  result.theColumns[0][0] = theColumns[0][0] + op.theColumns[0][0];
  result.theColumns[0][1] = theColumns[0][1] + op.theColumns[0][1];
  result.theColumns[1][0] = theColumns[1][0] + op.theColumns[1][0];
  result.theColumns[1][1] = theColumns[1][1] + op.theColumns[1][1];
  return result;
}

template<>
inline StaticMatrix<float, 2, 2> StaticMatrix<float, 2, 2>::operator*(const float& s) const {
  StaticMatrix<float, 2, 2> result;
  result.theColumns[0][0] = theColumns[0][0] * s;
  result.theColumns[0][1] = theColumns[0][1] * s;
  result.theColumns[1][0] = theColumns[1][0] * s;
  result.theColumns[1][1] = theColumns[1][1] * s;
  return result;
}

template<>
inline StaticVector<float, 4> StaticMatrix<float, 4, 4>::operator*(const StaticVector<float, 4>& v) const {
  StaticVector<float, 4> result;
  result[0] = theColumns[0][0] * v[0] + theColumns[1][0] * v[1] + theColumns[2][0] * v[2] + theColumns[3][0] * v[3];
  result[1] = theColumns[0][1] * v[0] + theColumns[1][1] * v[1] + theColumns[2][1] * v[2] + theColumns[3][1] * v[3];
  result[2] = theColumns[0][2] * v[0] + theColumns[1][2] * v[1] + theColumns[2][2] * v[2] + theColumns[3][2] * v[3];
  result[3] = theColumns[0][3] * v[0] + theColumns[1][3] * v[1] + theColumns[2][3] * v[2] + theColumns[3][3] * v[3];
  return result;
}
template<>
inline StaticVector<float, 3> StaticMatrix<float, 3, 3>::operator*(const StaticVector<float, 3>& v) const {
  StaticVector<float, 3> result;
  result[0] = theColumns[0][0] * v[0] + theColumns[1][0] * v[1] + theColumns[2][0] * v[2];
  result[1] = theColumns[0][1] * v[0] + theColumns[1][1] * v[1] + theColumns[2][1] * v[2];
  result[2] = theColumns[0][2] * v[0] + theColumns[1][2] * v[1] + theColumns[2][2] * v[2];
  return result;
}
template<>
inline StaticVector<float, 2> StaticMatrix<float, 2, 2>::operator*(const StaticVector<float, 2>& v) const {
  StaticVector<float, 2> result;
  result[0] = theColumns[0][0] * v[0] + theColumns[1][0] * v[1];
  result[1] = theColumns[0][1] * v[0] + theColumns[1][1] * v[1];
  return result;
}
template<>
inline StaticVector<float, 2> StaticMatrix<float, 3, 2>::operator*(const StaticVector<float, 3>& v) const {
  StaticVector<float, 2> result;
  result[0] = theColumns[0][0] * v[0] + theColumns[1][0] * v[1] + theColumns[2][0] * v[2];
  result[1] = theColumns[0][1] * v[0] + theColumns[1][1] * v[1] + theColumns[2][1] * v[2];
  return result;
}

template<>
inline float normQuad(StaticVector<float, 3> v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
};

template<>
inline float normQuad(StaticVector<float, 2> v) {
  return v[0] * v[0] + v[1] * v[1];
};


// --- partial specialization (efficiency)


template <class FloatTypeDest, class FloatTypeSource>
StaticVector<FloatTypeDest, 2> convertVector(const StaticVector<FloatTypeSource, 2>& source) {
  result[0] = (FloatTypeDest)source[0];
  result[1] = (FloatTypeDest)source[1];
  return result;
}

template <class FloatTypeDest, class FloatTypeSource>
StaticVector<FloatTypeDest, 3> convertVector(const StaticVector<FloatTypeSource, 3>& source) {
  result[0] = (FloatTypeDest)source[0];
  result[1] = (FloatTypeDest)source[1];
  result[2] = (FloatTypeDest)source[2];
  return result;
}

template <class FloatTypeDest, class FloatTypeSource>
StaticVector<FloatTypeDest, 4> convertVector(const StaticVector<FloatTypeSource, 4>& source) {
  result[0] = (FloatTypeDest)source[0];
  result[1] = (FloatTypeDest)source[1];
  result[2] = (FloatTypeDest)source[2];
  result[3] = (FloatTypeDest)source[3];
  return result;
}

template<>
inline float norm(StaticVector<float, 3> v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
};

template<>
inline float norm(StaticVector<float, 2> v) {
  return sqrt(v[0] * v[0] + v[1] * v[1]);
};

template<>
inline StaticVector<float, 1> makeVector<float, 1>(const float* data)
{
  StaticVector<float, 1> ret;
  ret[0] = data[0];
  return ret;
}


template<>
inline StaticVector<float, 2> makeVector<float, 2>(const float* data)
{
  StaticVector<float, 2> ret;
  ret[0] = data[0];
  ret[1] = data[1];
  return ret;
}

template<>
inline StaticVector<float, 3> makeVector<float, 3>(const float* data)
{
  StaticVector<float, 3> ret;
  ret[0] = data[0];
  ret[1] = data[1];
  ret[2] = data[2];
  return ret;
}

template<>
inline StaticVector<float, 4> makeVector<float, 4>(const float* data)
{
  StaticVector<float, 4> ret;
  ret[0] = data[0];
  ret[1] = data[1];
  ret[2] = data[2];
  ret[3] = data[3];
  return ret;
}


#endif  // end of explicit loop unrolling




// ----  helper routines


inline Vector2f makeVector2f(float x, float y) {
  Vector2f v;
  v[0] = x; v[1] = y;
  return v;
}

inline Vector2d makeVector2d(double x, double y) {
  Vector2d v;
  v[0] = x; v[1] = y;
  return v;
}

inline Vector3f makeVector3f(float x, float y, float z) {
  Vector3f v;
  v[0] = x; v[1] = y; v[2] = z;
  return v;
}

inline Vector3d makeVector3d(double x, double y, double z) {
  Vector3d v;
  v[0] = x; v[1] = y; v[2] = z;
  return v;
}

inline Vector4f makeVector4f(float x, float y, float z, float w) {
  Vector4f v;
  v[0] = x; v[1] = y; v[2] = z; v[3] = w;
  return v;
}

inline Vector4d makeVector4d(double x, double y, double z, double w) {
  Vector4d v;
  v[0] = x; v[1] = y; v[2] = z; v[3] = w;
  return v;
}

template<class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> makeVector(const FloatType* data) {
  StaticVector<FloatType, dim> ret;
  for (int i = 0; i < dim; i++)
    ret[i] = data[i];

  return ret;
}

inline Vector2f makeVector2fv(const float* data) {
  return makeVector<float, 2>(data);
}

inline Vector3f makeVector3fv(const float* data) {
  return makeVector<float, 3>(data);
}

inline Vector4f makeVector4fv(const float* data) {
  return makeVector<float, 4>(data);
}


inline Matrix4f makeMatrix4f(float m00, float m10, float m20, float m30,
  float m01, float m11, float m21, float m31,
  float m02, float m12, float m22, float m32,
  float m03, float m13, float m23, float m33) {
  Matrix4f m;
  m[0][0] = m00;  m[1][0] = m10;  m[2][0] = m20;  m[3][0] = m30;
  m[0][1] = m01;  m[1][1] = m11;  m[2][1] = m21;  m[3][1] = m31;
  m[0][2] = m02;  m[1][2] = m12;  m[2][2] = m22;  m[3][2] = m32;
  m[0][3] = m03;  m[1][3] = m13;  m[2][3] = m23;  m[3][3] = m33;
  return m;
}

inline Matrix3f makeMatrix3f(float m00, float m10, float m20,
  float m01, float m11, float m21,
  float m02, float m12, float m22) {
  Matrix3f m;
  m[0][0] = m00;  m[1][0] = m10;  m[2][0] = m20;
  m[0][1] = m01;  m[1][1] = m11;  m[2][1] = m21;
  m[0][2] = m02;  m[1][2] = m12;  m[2][2] = m22;
  return m;
}

inline Matrix3d makeMatrix3d(double m00, double m10, double m20,
  double m01, double m11, double m21,
  double m02, double m12, double m22) {
  Matrix3d m;
  m[0][0] = m00;  m[1][0] = m10;  m[2][0] = m20;
  m[0][1] = m01;  m[1][1] = m11;  m[2][1] = m21;
  m[0][2] = m02;  m[1][2] = m12;  m[2][2] = m22;
  return m;
}

inline Matrix2f makeMatrix2f(float m00, float m01,
  float m10, float m11) {
  Matrix2f m;
  m[0][0] = m00;  m[1][0] = m10;
  m[0][1] = m01;  m[1][1] = m11;
  return m;
}

inline Vector4f expand3To4(Vector3f v) {
  return makeVector4f(v[0], v[1], v[2], 1.0f);
}

inline Matrix4f expand3To4(Matrix3f v) {
  Matrix4f result;
  for (int x = 0; x < 3; x++) {
    for (int y = 0; y < 3; y++) {
      result[x][y] = v[x][y];
    }
  }
  for (int i = 0; i < 3; i++) {
    result[3][i] = 0.0;
  }
  for (int i = 0; i < 3; i++) {
    result[i][3] = 0.0;
  }
  result[3][3] = 1.0;
  return result;
}

inline Vector3f shrink4To3(Vector4f v) {
  return makeVector3f(v[0], v[1], v[2]);
}
inline Matrix3f shrink4To3(Matrix4f v) {
  Matrix3f result;
  for (int x = 0; x < 3; x++) {
    for (int y = 0; y < 3; y++) {
      result[x][y] = v[x][y];
    }
  }
  return result;
}

inline Vector3f projectHomogeneous4To3(Vector4f v) {
  return makeVector3f(v[0] / v[3], v[1] / v[3], v[2] / v[3]);
}

template <class FloatType, unsigned dim>
inline FloatType normQuad(StaticVector<FloatType, dim> v) {
  FloatType result = 0;
  for (int i = 0; i < dim; i++) {
    result += v[i] * v[i];
  }
  return result;
};

template <class FloatType, unsigned dim>
inline FloatType norm(StaticVector<FloatType, dim> v) {
  return sqrt(normQuad(v));
};

template <class FloatType, unsigned dim>
inline StaticVector<FloatType, dim> normalize(const StaticVector<FloatType, dim>& v) {
  return v / norm(v);
}

template<class FloatType, unsigned size>
inline unsigned searchPivot(const StaticMatrix<FloatType, size, size>& mat, unsigned step, bool* success, FloatType pivotEps) {
  unsigned result = step;
  FloatType max = fabs(mat[step][step]);
  for (unsigned p = step + 1; p < size; p++) {
    if (fabs(mat[step][p]) > max) {
      max = fabs(mat[step][p]);
      result = p;
    }
  }
  if (max <= pivotEps) {
    if (success != NULL) {
      *success = false;
      return 0;
    }
    else {
      
    }
  }
  else {
    if (success != NULL) {
      *success = true;
    }
    return result;
  }
};

template<class FloatType, unsigned size>
inline StaticMatrix<FloatType, size, size> invertMatrix(const StaticMatrix<FloatType, size, size>& mat, bool* success, FloatType pivotEps) {
  StaticMatrix<FloatType, size, size> result;
  StaticMatrix<FloatType, size, size> mCopy = mat;
  for (unsigned step = 0; step < size; step++) {
    unsigned pivot = searchPivot(mCopy, step, success, pivotEps);
    if (success != NULL) {
      if (!(*success)) {
        return result;
      }
    }
    if (step != pivot) {
      mCopy.changeRows(step, pivot);
      result.changeRows(step, pivot);
    }
//    if (mCopy[step][step] == 0) 
//        throw PException("invertMatrix: Matrix is singular.");
    FloatType m = (FloatType)1.0 / mCopy[step][step];
    mCopy.multRow(step, m);
    result.multRow(step, m);
    for (unsigned row = step + 1; row < size; row++) {
      FloatType m = -mCopy[step][row];
      mCopy.combineRows(row, step, m);
      result.combineRows(row, step, m);
    }
  }
  for (int step = size - 1; step >= 0; step--) {
    for (int row = step - 1; row >= 0; row--) {
      FloatType m = -mCopy[step][row];
      //mCopy.combineRows(row, step, m);
      result.combineRows(row, step, m);
    }
  }
  if (success != NULL) {
    *success = true;
  }
  return result;
};

inline void CompleteCoordinateFrame(const Vector3f& axis_z,
  Vector3f* axis_x,
  Vector3f* axis_y) {
  for (int i = 0; i < 3; i++) {
    (*axis_y)[i] = 0.0;
  }
  for (int i = 0; i < 3; i++) {
    if (fabs(axis_z[i]) < 0.6f) {
      (*axis_y)[i] = 1.f;
      break;
    }
  }
  (*axis_y) -= axis_z * (axis_z * (*axis_y));
  (*axis_y) /= float(sqrt(axis_y->getSqrNorm()));
  *axis_x = axis_y->crossProduct(axis_z);
}

inline float TriangleArea(const Vector3f& p1,
  const Vector3f& p2,
  const Vector3f& p3) {
  double l0 = (p1 - p2).getSqrNorm();
  double l1 = (p2 - p3).getSqrNorm();
  double l2 = (p3 - p1).getSqrNorm();

  return float(sqrt(2 * (l0 * l1 + l2 * (l0 + l1)) - l0 * l0 - l1 * l1 - l2 * l2) / 4);
}


inline double VectorAngle(const Vector3f& vec,
  const Vector3f& axis_x,
  const Vector3f& axis_y) {
  Vector2d direc;
  direc[0] = axis_x * vec;
  direc[1] = axis_y * vec;
  direc.normalize();

  double angle = acos(direc[0]);
  if (direc[1] < 0.f)
    angle = 4 * acos(0.0) - angle;
  return angle;
}

inline double VectorAngle(const Vector3d& vec,
  const Vector3d& axis_x,
  const Vector3d& axis_y) {
  Vector2d direc;
  direc[0] = axis_x * vec;
  direc[1] = axis_y * vec;
  direc.normalize();

  double angle = acos(direc[0]);
  if (direc[1] <= 0)
    angle = 4 * acos(0.0) - angle;
  return angle;
}

inline double AngleBetweenVectors(const Vector3f& origin,
  const Vector3f& v1,
  const Vector3f& v2) {
  double a = (v1 - origin).getSqrNorm();
  double b = (v2 - origin).getSqrNorm();
  double c = (v1 - v2).getSqrNorm();

  double cos_angle = (a + b - c) / 2 / sqrt(a * b);
  return acos(cos_angle);
}

inline double AngleBetweenVectors(const Vector3d& origin,
  const Vector3d& v1,
  const Vector3d& v2) {
  double a = (v1 - origin).getSqrNorm();
  double b = (v2 - origin).getSqrNorm();
  double c = (v1 - v2).getSqrNorm();

  double cos_angle = (a + b - c) / 2 / sqrt(a * b);
  return acos(cos_angle);
}

#endif
