#ifndef simul_reg_and_recons_h_
#define simul_reg_and_recons_h_

#include "linear_algebra.h"
#include "point_cloud.h"
#include "octree.h"

#include <vector>
using namespace std;


struct SimulRegAndReconsPara {
public:
  SimulRegAndReconsPara() {
    gridSize = 0.01;
    stride = 1;
    weightPoint2PlaneDis = 0.9;
    minNumPointsPerCell = 4;
    maxNumPointsPerCell = 256;
    numAlternatingIterations = 6;
  }
  ~SimulRegAndReconsPara() {
  }
  double gridSize; // The size of the grid we use to build the reference surface for alignment
  int stride; //For effiency, we only align points in cells specified by stride 'stride'
  double weightPoint2PlaneDis; // weight of the point-plane distance
  int minNumPointsPerCell; // Minimum number of points per cell
  int maxNumPointsPerCell; // Maximum number of points per cell used to perform PCA analysis
  int numAlternatingIterations; //Number of alternating iterations
};

struct PointCorres {
  PointCorres() {
    sourcePointId = 0;
    targetPointId = 0;
    weight = 0.f;
  }
  ~PointCorres() {}
  unsigned sourcePointId;
  unsigned targetPointId;
  float weight;
};

// Compute the reference surface from the input scans
// Find the target point of each scan point
class SimulRegAndRecons {
 public:
   SimulRegAndRecons() {
     max_depth_ = 0;
     depth_ = 0;
     cell_offsets_[0] = cell_offsets_[1] = cell_offsets_[2] = 0;
  }
  ~SimulRegAndRecons() {
  }
  void AlternatingOpt(const vector<PointCloud> &scans,
    const SimulRegAndReconsPara &para,
    vector<Affine3d> *opt_poses);
  void GenerateSuperScan(const vector<PointCloud> &scans,
    const vector<Affine3d> &opt_poses, 
    const SimulRegAndReconsPara &para,
    PointCloud *super_scan);
  void LatentSurfOpt(const vector<PointCloud>& scans,
    const vector<Affine3d>& opt_poses,
    const SimulRegAndReconsPara& para,
    vector<vector<PointCorres>>* corres,
    PointCloud* latentSurf);
 private:
   Octree3D* GenerateOctree(
     const vector<Surfel3D>& points,
     const double& grid_size,
     const unsigned& stride);
   void InsertAPoint(const int& point_index,
     Node3D** current_node);
   void SurfelFitting(const vector<Vector3f>& poss,
     const vector<Vector3f>& nors,
     const float& weight_nor,
     Surfel3D* fit_sur);
 private:
   int max_depth_;
   int depth_;
   int cell_offsets_[3];
};

#endif