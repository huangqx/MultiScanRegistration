//---------------------------------------------------------------------------
#ifndef linear_algebra_h_
#define linear_algebra_h_

class Vector3f {
 private:
   float v[3];
 public:
   float& operator[](const unsigned &index);
   const float& operator[](const unsigned &index) const;

   //  ----  operators +, -, *, /, crossProduct
   Vector3f operator+(const Vector3f &op) const;
   Vector3f operator-(const Vector3f &op) const;
   Vector3f operator-() const;
   float operator*(const Vector3f &op) const;
   Vector3f operator*(const float &s) const;
   Vector3f operator/(const float &s) const;
   Vector3f crossProduct(const Vector3f &op) const;
   Vector3f componentProduct(const Vector3f &op) const;

   //  ---- operators +=, -=, *=, /=
   Vector3f operator+=(const Vector3f &op);
   Vector3f operator-=(const Vector3f &op);
   Vector3f operator*=(const float &s);
   Vector3f operator/=(const float &s);


   void normalize();
   float getSqrNorm() const;
};

class Vector3d {
private:
  double v[3];
public:
  double& operator[](const unsigned &index);
  const double& operator[](const unsigned &index) const;

  //  ----  operators +, -, *, /, crossProduct
  Vector3d operator+(const Vector3d &op) const;
  Vector3d operator-(const Vector3d &op) const;
  Vector3d operator-() const;
  double operator*(const Vector3d &op) const;
  Vector3d operator*(const double &s) const;
  Vector3d operator/(const double &s) const;
  Vector3d crossProduct(const Vector3d &op) const;
  Vector3d componentProduct(const Vector3d &op) const;

  //  ---- operators +=, -=, *=, /=
  Vector3d operator+=(const Vector3d &op);
  Vector3d operator-=(const Vector3d &op);
  Vector3d operator*=(const double &s);
  Vector3d operator/=(const double &s);


  void normalize();
  double getSqrNorm() const;
};

#endif
