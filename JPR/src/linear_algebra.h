//---------------------------------------------------------------------------
#ifndef linear_algebra_h_
#define linear_algebra_h_

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


#endif
