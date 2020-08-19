//---------------------------------------------------------------------------
#include "linear_algebra.h"
#include "math.h"
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------



float& Vector3f::operator[](const unsigned &index) {
  return v[index];
}

const float& Vector3f::operator[](const unsigned &index) const{
  return v[index];
}

Vector3f Vector3f::operator+(const Vector3f &op) const {
  Vector3f result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] + op.v[i];
  }
  return result;
}

Vector3f Vector3f::operator-(const Vector3f &op) const {
  Vector3f result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] - op.v[i];
  }
  return result;
}

Vector3f Vector3f::operator-() const {
  Vector3f result;
  for (unsigned i=0; i<3; i++) {
    result[i] = -v[i];
  }
  return result;
}

float Vector3f::operator*(const Vector3f &op) const {
  float result = 0.0;
  for (unsigned i=0; i<3; i++) {
    result += v[i] * op.v[i];
  }
  return result;
}

Vector3f Vector3f::operator*(const float &s) const {
  Vector3f result;
  for (unsigned i=0; i<3; i++) {
    result[i] = s * v[i];
  }
  return result;
}

Vector3f Vector3f::operator/(const float &s) const {
  Vector3f result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] / s;
  }
  return result;
}

Vector3f Vector3f::crossProduct(const Vector3f &op) const {
  Vector3f result;
  result[0] = v[1]*op.v[2] - v[2]*op.v[1];
  result[1] = v[2]*op.v[0] - v[0]*op.v[2];
  result[2] = v[0]*op.v[1] - v[1]*op.v[0];
  return result;
}

Vector3f Vector3f::componentProduct(const Vector3f &op) const {
  Vector3f result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] * op.v[i];
  }
  return result;
}

Vector3f Vector3f::operator+=(const Vector3f &op) {
  for (unsigned i=0; i<3; i++) {
    v[i] += op.v[i];
  }
  return *this;
}

Vector3f Vector3f::operator-=(const Vector3f &op) {
  for (unsigned i=0; i<3; i++) {
    v[i] -= op.v[i];
  }
  return *this;
}

Vector3f Vector3f::operator*=(const float &s) {
  for (unsigned i=0; i<3; i++) {
    v[i] *= s;
  }
  return *this;
}

Vector3f Vector3f::operator/=(const float &s) {
  for (unsigned i=0; i<3; i++) {
    v[i] /= s;
  }
  return *this;
}

void Vector3f::normalize() {
  float norm = 0.0;
  for (int i=0; i<3; i++) {
    norm += v[i] * v[i];
  }
  if ( norm > 1e-12f )
    (*this) = (*this) / sqrt( norm );
}

float Vector3f::getSqrNorm() const {
  float norm2 = 0.f;
  for (int i=0; i<3; i++) {
    norm2 += v[i] * v[i];
  }
  return norm2;
}


/************************************************************************/
/*                                                                      */
/************************************************************************/

double& Vector3d::operator[](const unsigned &index) {
  return v[index];
}

const double& Vector3d::operator[](const unsigned &index) const{
  return v[index];
}

Vector3d Vector3d::operator+(const Vector3d &op) const {
  Vector3d result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] + op.v[i];
  }
  return result;
}

Vector3d Vector3d::operator-(const Vector3d &op) const {
  Vector3d result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] - op.v[i];
  }
  return result;
}

Vector3d Vector3d::operator-() const {
  Vector3d result;
  for (unsigned i=0; i<3; i++) {
    result[i] = -v[i];
  }
  return result;
}

double Vector3d::operator*(const Vector3d &op) const {
  double result = 0.0;
  for (unsigned i=0; i<3; i++) {
    result += v[i] * op.v[i];
  }
  return result;
}

Vector3d Vector3d::operator*(const double &s) const {
  Vector3d result;
  for (unsigned i=0; i<3; i++) {
    result[i] = s * v[i];
  }
  return result;
}

Vector3d Vector3d::operator/(const double &s) const {
  Vector3d result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] / s;
  }
  return result;
}

Vector3d Vector3d::crossProduct(const Vector3d &op) const {
  Vector3d result;
  result[0] = v[1]*op.v[2] - v[2]*op.v[1];
  result[1] = v[2]*op.v[0] - v[0]*op.v[2];
  result[2] = v[0]*op.v[1] - v[1]*op.v[0];
  return result;
}

Vector3d Vector3d::componentProduct(const Vector3d &op) const {
  Vector3d result;
  for (unsigned i=0; i<3; i++) {
    result[i] = v[i] * op.v[i];
  }
  return result;
}

Vector3d Vector3d::operator+=(const Vector3d &op) {
  for (unsigned i=0; i<3; i++) {
    v[i] += op.v[i];
  }
  return *this;
}

Vector3d Vector3d::operator-=(const Vector3d &op) {
  for (unsigned i=0; i<3; i++) {
    v[i] -= op.v[i];
  }
  return *this;
}

Vector3d Vector3d::operator*=(const double &s) {
  for (unsigned i=0; i<3; i++) {
    v[i] *= s;
  }
  return *this;
}

Vector3d Vector3d::operator/=(const double &s) {
  for (unsigned i=0; i<3; i++) {
    v[i] /= s;
  }
  return *this;
}

void Vector3d::normalize() {
  double norm = 0.0;
  for (int i=0; i<3; i++) {
    norm += v[i] * v[i];
  }
  if ( norm > 1e-16 )
    (*this) = (*this) / sqrt( norm );
}

double Vector3d::getSqrNorm() const {
  double norm2 = 0.f;
  for (int i=0; i<3; i++) {
    norm2 += v[i] * v[i];
  }
  return norm2;
}
