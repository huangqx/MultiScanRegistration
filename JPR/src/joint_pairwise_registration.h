#ifndef joint_pairwise_registration_h_
#define joint_pairwise_registration_h_

#include "linear_algebra.h"
#include "dynamic_linear_algebra.h"
#include "point_cloud.h"
#include "octree.h"
#include "ANN.h"

#include <vector>
using namespace std;


class ANNDataStructure {
 public:
   ANNDataStructure() {
     dataPts_fixed_ = NULL;
     kdTree_fixed_ = NULL;
     queryPt_ = NULL;
     nnIdx_ = NULL;
     dists_ = NULL;
   }
   ~ANNDataStructure() {
     ClearANN();
   }
   void InitializeANN(const PointCloud& pc,
     const double& weightNormal,
     const double& weightColor);
   void ClearANN();
   ANNpointArray	dataPts_fixed_;
   ANNkd_tree* kdTree_fixed_;
   ANNpoint queryPt_;
   ANNidxArray	nnIdx_;
   ANNdistArray dists_;
};

// Compute the reference surface from the input scans
// Find the target point of each scan point
class JointPairwiseRegistration {
 public:
  JointPairwiseRegistration() {
  }
  ~JointPairwiseRegistration() {
  }
  void Compute(const vector<PointCloud> &scans,
    const float& max_distance,
    const float &overlap_ratio,
    const float &weight_normal,
    const float &weight_color,
    const unsigned &down_sampling_rate,
    const unsigned &num_reweighting_iters,
    const unsigned &num_gauss_newton_iters, 
    const float &weight_point2planeDis,
    vector<Affine3d> *opt_poses,
    vector<PointCorres2> *pointcorres);
 private:
  void Registration_With_KnownCorres(const vector<PointCloud>& scans,
    const vector<PointCorres2>& pointcorres,
    const unsigned& num_gauss_newton_iters,
    const float& weight_point2planeDis,
    vector<Affine3d>* opt_poses);
  double RegistrationObjectiveValue(const vector<PointCloud>& scans,
    const vector<PointCorres2>& pointcorres,
    const vector<Affine3d>& opt_poses,
    const float& weight_point2planeDis);
  void Initialize_ANN(const vector<PointCloud>& scans,
    const float &weightNormal,
    const float &weightColor);
  void Clear_ANN();
  void Correspondence_Reweighting(vector<PointCorres2>* pointcorres);
  void Compute_Correspondences(const vector<PointCloud>& scans,
    const vector<Affine3d>& cur_poses,
    const float& max_distance,
    const float& overlap_ratio,
    const float& weight_normal,
    const float& weight_color,
    const unsigned& down_sampling_rate,
    vector<PointCorres2>* pointcorres);
  bool SolveNxN(const DMatrixD& matA, const DVectorD& vecb,
    DVectorD* vecx);
  vector<ANNDataStructure> nn_search_dss_;
};

#endif