#include "affine_transformation.h"

#include <math.h>
#include <string>

using namespace std;

inline Vector2d& Affine2d::operator[](const unsigned& index) {
  return v[index];
}
inline const Vector2d& Affine2d::operator[](const unsigned& index) const {
  return v[index];
}

inline Affine2d Affine2d::operator+(const Affine2d& op) const {
  Affine2d result;
  result.v[0] = v[0] + op.v[0];
  result.v[1] = v[1] + op.v[1];
  result.v[2] = v[2] + op.v[2];
  return result;
}


inline Affine2d Affine2d::operator-(const Affine2d& op) const {
  Affine2d result;
  result.v[0] = v[0] - op.v[0];
  result.v[1] = v[1] - op.v[1];
  result.v[2] = v[2] - op.v[2];
  return result;
}

inline Affine2d Affine2d::operator-() const {
  Affine2d result;
  result.v[0] = -v[0];
  result.v[1] = -v[1];
  result.v[2] = -v[2];
  return result;
}

inline Affine2d Affine2d::operator*(const double& s) const {
  Affine2d result;
  result.v[0] = v[0] * s;
  result.v[1] = v[1] * s;
  result.v[2] = v[2] * s;
  return result;
}

inline Affine2d Affine2d::operator/(const double& s) const {
  Affine2d result;
  result.v[0] = v[0] / s;
  result.v[1] = v[1] / s;
  result.v[2] = v[2] / s;
  return result;
}

inline Affine2d Affine2d::operator+=(const Affine2d& op) {
  v[0] += op.v[0];
  v[1] += op.v[1];
  v[2] += op.v[2];
  return *this;
}

inline Affine2d Affine2d::operator-=(const Affine2d& op) {
  v[0] -= op.v[0];
  v[1] -= op.v[1];
  v[2] -= op.v[2];
  return *this;
}

inline Affine2d Affine2d::operator*=(const double& op) {
  v[0] *= op;
  v[1] *= op;
  v[2] *= op;
  return *this;
}

inline Affine2d Affine2d::operator/=(const double& op) {
  v[0] /= op;
  v[1] /= op;
  v[2] /= op;
  return *this;
}

inline Affine2d Affine2d::operator=(const Affine2d& op) {
  v[0] = op.v[0];
  v[1] = op.v[1];
  v[2] = op.v[2];
  return *this;
}

inline Vector2d Affine2d::operator*(const Vector2d& vec) const {
  return v[0] + v[1] * vec[0] + v[2] * vec[1];
}

inline Affine2d::Affine2d() {
  Initialize();
}

inline void Affine2d::Initialize() {
  SetZero();
  v[1][0] = 1.0;
  v[2][1] = 1.0;
}

inline void Affine2d::SetZero() {
  v[0][0] = v[0][1] = 0;
  v[1][0] = v[1][1] = 0;
  v[2][0] = v[2][1] = 0;
}

inline double Affine2d::Det() {
  return v[1][0] * v[2][1] - v[2][0] * v[1][1];
}

