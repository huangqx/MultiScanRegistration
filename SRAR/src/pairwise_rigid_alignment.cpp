#include "pairwise_rigid_alignment.h"
#include "dynamic_linear_algebra.h"
#include <algorithm>

void PairwiseRigidAlign::Compute_KnownCorres(const PointCloud& fixedSurf,
  const PointCloud& movingSurf,
  const vector<int>& movingPointIds,
  const vector<int>& fixedPointIds,
  const vector<double>& corresWeights,
  const double& weightPoint2Plane,
  const int & maxIters,
  Affine3d* rigidTransform) {
  // Precompute the Hessian matrix in the affine transformation space
  Matrix12d hessian_aff_pn;
  Vector12d gradient_aff_pn;
  Matrix4d hessian_aff_pp;
  Vector12d gradient_aff_pp;
  double J_pp[4];
  double J_pn[12];
  double g_pp[3];
  double g_pn;
  for (int i = 0; i < 12; ++i) {
    for (int j = i; j < 12; ++j)
      hessian_aff_pn[i][j] = 0.0;
    gradient_aff_pn[i] = 0.0;
    gradient_aff_pp[i] = 0.0;
    J_pn[i] = 0.0;
  }
  for (int i = 0; i < 4; ++i) {
    for (int j = i; j < 4; ++j)
      hessian_aff_pp[i][j] = 0.0;
    J_pp[i] = 0.0;
  }

  Vector3d footPos, // Target position of each moving point
    footNor; // Normal at the target position
  for (unsigned i = 0; i < movingPointIds.size(); ++i) {
    double weight = corresWeights[i];
    const Surfel3D& point = (*movingSurf.GetPointArray())[movingPointIds[i]];
    
    const Surfel3D& foot = (*fixedSurf.GetPointArray())[fixedPointIds[i]];
    for (int i = 0; i < 3; ++i) {
      footPos[i] = foot.position[i];
      footNor[i] = foot.normal[i];
    }
    // Update the term related to point-2-plane distance
    g_pn = footPos * footNor;
    for (int i = 0; i < 3; ++i) {
      J_pn[i] = footNor[i];
      for (int j = 0; j < 3; ++j) {
        J_pn[i+3*j+3] = footNor[i]*point.position[j];
      }
    }
    for (int i = 0; i < 12; ++i) {
      gradient_aff_pn[i] += J_pn[i] * g_pn * weight;
      for (int j = i; j < 12; ++j)
        hessian_aff_pn[i][j] += J_pn[i] * J_pn[j] * weight;
    }
    // Update the term related to point-2-point distance
    J_pp[0] = 1;
    J_pp[1] = point.position[0];
    J_pp[2] = point.position[1];
    J_pp[3] = point.position[2];
    g_pp[0] = footPos[0];
    g_pp[1] = footPos[1];
    g_pp[2] = footPos[2];
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 3; ++j)
        gradient_aff_pp[3*i+j] += J_pp[i]*g_pp[j]* weight;
      for (int j = i; j < 4; ++j)
        hessian_aff_pp[i][j] += J_pp[i] * J_pp[j]* weight;
    }
  }
  Matrix12d hessian_aff = hessian_aff_pn*weightPoint2Plane;
  double weightPoint2Point = 1 - weightPoint2Plane;
  for (int i = 0; i < 4; ++i)
    for (int j = i; j < 4; ++j)
      for (int k = 0; k < 3; ++k)
        hessian_aff[3*i+k][3*j+k] += hessian_aff_pp[i][j] * weightPoint2Point;
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < i; ++j)
      hessian_aff[i][j] = hessian_aff[j][i];

  Vector12d gradient_aff = gradient_aff_pn*weightPoint2Plane 
    + gradient_aff_pp*weightPoint2Point;

  for (int iter = 0; iter < maxIters; ++iter) {
    Vector12d cur_sol;
    rigidTransform->GetData(&cur_sol);

    double J[9][3];
    for (int i = 0; i < 3; ++i) {
      Vector3d vec = (*rigidTransform)[i+1];
      J[3*i][0] = 0.0;
      J[3*i][1] = vec[2];
      J[3*i][2] = -vec[1];
      J[3*i+1][0] = -vec[2];
      J[3*i+1][1] = 0.0;
      J[3*i+1][2] = vec[0];
      J[3*i+2][0] = vec[1];
      J[3*i+2][1] = -vec[0];
      J[3*i+2][2] = 0.0;
    }

    Vector12d g = gradient_aff - hessian_aff * cur_sol;
    Matrix6d hessian;
    Vector6d gradient;
    for (int i = 0; i < 3; ++i) {
      // hessian
      for (int j = i; j < 3; ++j) {
        hessian[i][j] = hessian_aff[i][j];
        int i1 = i + 3;
        int j1 = j + 3;
        hessian[i1][j1] = 0.0;
        for (int k = 0; k < 9; ++k)
          for (int k1 = 0; k1 < 9; ++k1)
            hessian[i1][j1] += J[k][i] * J[k1][j] * hessian_aff[k + 3][k1 + 3];
      }
      for (int j = 0; j < 3; ++j) {
        int j1 = j + 3;
        hessian[i][j1] = 0.0;
        for (int k = 0; k < 9; ++k)
          hessian[i][j1] += J[k][j] * hessian_aff[i][k + 3];
      }
      // gradient
      gradient[i] = g[i];
      gradient[i + 3] = 0.0;
      for (int k = 0; k < 9; ++k)
        gradient[i + 3] += g[3 + k] * J[k][i];
    }
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < i; ++j)
        hessian[i][j] = hessian[j][i];

    double velocity[6];
    Solve6x6(hessian, gradient, velocity);
    double norm2 = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]
      + velocity[3] * velocity[3] + velocity[4] * velocity[4] + velocity[5] * velocity[5]);
    if (norm2 < 1e-6)
      break;
    Affine3d newMotion(velocity);
    Affine3d rigid_next = newMotion*(*rigidTransform);
    Vector12d next_sol;
    rigid_next.GetData(&next_sol);
    double e_cur = cur_sol * ((hessian_aff * cur_sol)*0.5 - gradient_aff);
    double e_next = next_sol * ((hessian_aff * next_sol) * 0.5 - gradient_aff);
    if (e_next < e_cur) {
      *rigidTransform = rigid_next;
    } else {
      double lambda = 0;
      for (int i = 0; i < 6; ++i)
        lambda += hessian[i][i]/6;
      lambda = lambda / 1e4;
      bool success = false;
      for (int linesearchId = 0; linesearchId < 15; linesearchId++) {
        for (int i = 0; i < 6; ++i)
          hessian[i][i] += lambda;
        Solve6x6(hessian, gradient, velocity);
        Affine3d newMotion(velocity);
        Affine3d rigid_next = newMotion * (*rigidTransform);
        Vector12d next_sol;
        rigid_next.GetData(&next_sol);
        double e_next = next_sol * ((hessian_aff * next_sol) * 0.5 - gradient_aff);
        if (e_next < e_cur) {
          *rigidTransform = rigid_next;
          success = true;
          break;
        }
        lambda = lambda * 2;
      }
      if (!success)
        break;
    }
  }
}


bool PairwiseRigidAlign::Solve6x6(const Matrix6d &hessian, Vector6d &gradient,
  double *velocity) {
  // Using LLT factorization to solve the symmetric linear system
  const double fTolerance = 1e-20;
  double afV[6], Lower[6][6];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j)
      Lower[i][j] = hessian[i][j];
    velocity[i] = gradient[i];
  }

  for (int i1 = 0; i1 < 6; ++i1) {
    for (int i0 = 0; i0 < i1; ++i0)
      afV[i0] = Lower[i1][i0]*Lower[i0][i0];

    afV[i1] = Lower[i1][i1];
    for (int i0 = 0; i0 < i1; ++i0)
      afV[i1] -= Lower[i1][i0]*afV[i0];

    Lower[i1][i1] = afV[i1];
    if (fabs(afV[i1]) <= fTolerance) //singular
      return false;

    double fInv = 1.0f/afV[i1];
    for (int i0 = i1+1; i0 < 6; ++i0) {
      for (int i2 = 0; i2 < i1; ++i2)
        Lower[i0][i1] -= Lower[i0][i2]*afV[i2];
      Lower[i0][i1] *= fInv;
    }
  }

  // Solve Ax = B.
  // Forward substitution
  for (int i0 = 0; i0 < 6; ++i0) {
    for (int i1 = 0; i1 < i0; ++i1)
      velocity[i0] -= Lower[i0][i1]*velocity[i1];
  }

  // Diagonal division:  Let y = L^t x, then Dy = z.  Algorithm stores
  // y terms in B vector.
  for (int i0 = 0; i0 < 6; ++i0) {
    if (fabs(Lower[i0][i0]) <= fTolerance )
      return false;
    velocity[i0] /= Lower[i0][i0];
  }

  // Back substitution:  Solve L^t x = y.  Algorithm stores x terms in
  // B vector.
  for (int i0 = 4; i0 >= 0; i0--) {
    for (int i1 = i0+1; i1 < 6; ++i1)
      velocity[i0] -= Lower[i1][i0]*velocity[i1];
  }
  return true;
}