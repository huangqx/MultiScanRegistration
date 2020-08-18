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

/************************************************************************/
/* Affine3d                                                             */
/************************************************************************/

inline Affine3d::Affine3d(double* velocity) {
  Vector3d barC, C;
  barC[0] = velocity[0];
  barC[1] = velocity[1];
  barC[2] = velocity[2];
  C[0] = velocity[3];
  C[1] = velocity[4];
  C[2] = velocity[5];
  double theta = sqrt(C * C);
  if (theta < 1e-8) {
    v[0] = barC;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        v[i + 1][j] = 0.0;
      }
      v[i + 1][i] = 1.0;
    }
  }
  else {
    C /= theta;
    double sintheta = sin(theta);
    double sintheta2 = sin(theta / 2);
    double c1 = sintheta2 * sintheta2 * 2 / theta;
    double c2 = 1 - (sintheta / theta);
    double c0 = 2 * sintheta2 * sintheta2;

    v[1][0] = (1 - c0) + c0 * C[0] * C[0];
    v[1][1] = sintheta * C[2] + c0 * C[0] * C[1];
    v[1][2] = -sintheta * C[1] + c0 * C[0] * C[2];

    v[2][0] = -sintheta * C[2] + c0 * C[1] * C[0];
    v[2][1] = (1 - c0) + c0 * C[1] * C[1];
    v[2][2] = sintheta * C[0] + c0 * C[1] * C[2];

    v[3][0] = sintheta * C[1] + c0 * C[2] * C[0];
    v[3][1] = -sintheta * C[0] + c0 * C[2] * C[1];
    v[3][2] = (1 - c0) + c0 * C[2] * C[2];

    Matrix3d matV;
    matV[0][0] = (1 - c2) + c2 * C[0] * C[0];
    matV[0][1] = c1 * C[2] + c2 * C[0] * C[1];
    matV[0][2] = -c1 * C[1] + c2 * C[0] * C[2];

    matV[1][0] = -c1 * C[2] + c2 * C[1] * C[0];
    matV[1][1] = (1 - c2) + c2 * C[1] * C[1];
    matV[1][2] = c1 * C[0] + c2 * C[1] * C[2];

    matV[2][0] = c1 * C[1] + c2 * C[2] * C[0];
    matV[2][1] = -c1 * C[0] + c2 * C[2] * C[1];
    matV[2][2] = (1 - c2) + c2 * C[2] * C[2];

    v[0] = matV * barC;
  }
}

inline Vector3d& Affine3d::operator[](const unsigned& index) {
  return v[index];
}

inline const Vector3d& Affine3d::operator[](const unsigned& index) const {
  return v[index];
}

inline Affine3d Affine3d::operator+(const Affine3d& op) const {
  Affine3d result;
  result.v[0] = v[0] + op.v[0];
  result.v[1] = v[1] + op.v[1];
  result.v[2] = v[2] + op.v[2];
  result.v[3] = v[3] + op.v[3];
  return result;
}

inline Affine3d Affine3d::operator-(const Affine3d& op) const {
  Affine3d result;
  result.v[0] = v[0] - op.v[0];
  result.v[1] = v[1] - op.v[1];
  result.v[2] = v[2] - op.v[2];
  result.v[3] = v[3] - op.v[3];
  return result;
}

inline Affine3d Affine3d::operator-() const {
  Affine3d result;
  result.v[0] = -v[0];
  result.v[1] = -v[1];
  result.v[2] = -v[2];
  result.v[3] = -v[3];
  return result;
}

inline Affine3d Affine3d::operator*(const double& s) const {
  Affine3d result;
  result.v[0] = v[0] * s;
  result.v[1] = v[1] * s;
  result.v[2] = v[2] * s;
  result.v[3] = v[3] * s;
  return result;
}

inline Affine3d Affine3d::operator*(const Affine3d& op) const {
  Affine3d result;
  result.v[0] = v[0] + v[1] * op.v[0][0] + v[2] * op.v[0][1] + v[3] * op.v[0][2];
  for (int k = 1; k < 4; k++)
    result.v[k] = v[1] * op.v[k][0] + v[2] * op.v[k][1] + v[3] * op.v[k][2];
  return result;
}

inline Affine3d Affine3d::operator/(const double& s) const {
  Affine3d result;
  result.v[0] = v[0] / s;
  result.v[1] = v[1] / s;
  result.v[2] = v[2] / s;
  result.v[3] = v[3] / s;
  return result;
}

inline Affine3d Affine3d::operator+=(const Affine3d& op) {
  v[0] += op.v[0];
  v[1] += op.v[1];
  v[2] += op.v[2];
  v[3] += op.v[3];
  return *this;
}

inline Affine3d Affine3d::operator-=(const Affine3d& op) {
  v[0] -= op.v[0];
  v[1] -= op.v[1];
  v[2] -= op.v[2];
  v[3] -= op.v[3];
  return *this;
}

inline Affine3d Affine3d::operator*=(const double& op) {
  v[0] *= op;
  v[1] *= op;
  v[2] *= op;
  v[3] *= op;
  return *this;
}

inline Affine3d Affine3d::operator/=(const double& op) {
  v[0] /= op;
  v[1] /= op;
  v[2] /= op;
  v[3] /= op;
  return *this;
}

inline Affine3d Affine3d::operator=(const Affine3d& op) {
  v[0] = op.v[0];
  v[1] = op.v[1];
  v[2] = op.v[2];
  v[3] = op.v[3];
  return *this;
}

inline Vector3d Affine3d::operator*(const Vector3d& vec) const {
  return v[0] + v[1] * vec[0] + v[2] * vec[1] + v[3] * vec[2];
}

inline Affine3d::Affine3d() {
  Initialize();
}

inline void Affine3d::Initialize() {
  SetZero();
  v[1][0] = 1.0;
  v[2][1] = 1.0;
  v[3][2] = 1.0;
}


inline void Affine3d::SetZero() {
  v[0][0] = v[0][1] = v[0][2] = 0;
  v[1][0] = v[1][1] = v[1][2] = 0;
  v[2][0] = v[2][1] = v[2][2] = 0;
  v[3][0] = v[3][1] = v[3][2] = 0;
}

inline double Affine3d::Det() {
  return v[1][0] * (v[2][1] * v[3][2] - v[3][1] * v[2][2])
    + v[2][0] * (v[3][1] * v[1][2] - v[1][1] * v[3][2])
    + v[3][0] * (v[1][1] * v[2][2] - v[2][1] * v[1][2]);
}

inline Affine3d	Affine3d::Inverse() {
  Affine3d	result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      result.v[i + 1][j] = v[j + 1][i];
    }
  }
  result.v[0] = -(result.v[1] * v[0][0]
    + result.v[2] * v[0][1]
    + result.v[3] * v[0][2]);
  return result;
}

inline void Affine3d::Normalize() {
  v[2] = v[2] - v[1] * (v[2] * v[1]);
  v[2].normalize();
  v[3] = v[3] - v[1] * (v[3] * v[1]) - v[2] * (v[3] * v[2]);
  v[3].normalize();
}

inline Affine3d	Affine3d::Inverse() const {
  Affine3d	result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      result.v[i + 1][j] = v[j + 1][i];
    }
  }
  result.v[0] = -(result.v[1] * v[0][0]
    + result.v[2] * v[0][1]
    + result.v[3] * v[0][2]);
  return result;
}

inline void Affine3d::GetData(Vector12d* vec) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
      (*vec)[3 * i + j] = v[i][j];
}

inline void Affine3d::SetData(const Vector12d& vec) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
      v[i][j] = vec[3 * i + j];
}
