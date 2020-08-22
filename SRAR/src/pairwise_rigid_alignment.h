#ifndef pairwise_rigid_alignment_h_
#define pairwise_rigid_alignment_h_

#include "point_cloud.h"
#include "affine_transformation.h"

#include <vector>
using namespace std;


struct PairMatch {
public:
  PairMatch() {
    fixedSurfId = 0;
    movingSurfId = 0;
    relativePose.Initialize();
    alignmentError = 0.f;
  }
  ~PairMatch() {

  }
  unsigned fixedSurfId;
  unsigned movingSurfId;
  Affine3d relativePose;
  float alignmentError;
};

struct PairwiseRigidAligPara {
 public:
  PairwiseRigidAligPara() {
    numIterations_reweighting = 6;
    numIterations_gauss_newton = 6;
    numSamples = 1024;
    weightPoint2PlaneDis = 0.95;
    weightNormal = 1.0;
    weightColor = 0.2;
    overlap_ratio = 0.25;
  }
  ~PairwiseRigidAligPara() {
  }
  double weightPoint2PlaneDis;
  double weightNormal;
  double weightColor;
  int numSamples;

  // For each level, we use a fixed estimate of the averaged
  // distance between points and their target points
  int numIterations_reweighting; // Number of reweighting steps
  int numIterations_gauss_newton; // Number of Gaussian-Newton steps

  double overlap_ratio; // used in determine the overlapping ratio
};


class PairwiseRigidAlign {
 public:
  PairwiseRigidAlign() {
      
  }
  ~PairwiseRigidAlign() {
      
  }
  void Compute_KnownCorres(const PointCloud& fixedSurf,
    const PointCloud& movingSurf,
    const vector<int>& movingPointIds,
    const vector<int>& fixedPointIds,
    const vector<double>& corresWeights,
    const double& weightPoint2Plane,
    const int& maxIters,
    Affine3d* rigidTransform);
 private:
  bool Solve6x6(const Matrix6d &hessian, Vector6d &gradient,
    double *velocity);
};

#endif