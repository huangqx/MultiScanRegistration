/*---
function [scanposes, corres] = jpr_interface(datapoints, scan_offsets, paras)
Input: 
Output
--*/

#include "mex.h"
#include "point_cloud.h"
#include "linear_algebra.h"
#include "affine_transformation.h"
#include "joint_pairwise_registration.h"
#include <vector>
#include <algorithm>
using namespace std;

// The three input matrices are
// 1) The mesh vertices
// 2) The mesh faces
// 3) The camera parameters

// The output parameters
// 1) The intersecting pixels

void mexFunction(
     int nargout,
     mxArray *output[],
     int nargin,
     const mxArray *input[]) {
  /* check argument */
  if (nargin != 3) {
    mexErrMsgTxt("Three input arguments required.");
  }
  if (nargout != 2) {
    mexErrMsgTxt("Incorrect number of output arguments."); 
  }

  double *pointData = (double*)mxGetData(input[0]);
  unsigned numPoints = static_cast<unsigned> (mxGetN(input[0]));
  double *offsetData = (double*)mxGetData(input[1]);
  unsigned numScans = static_cast<unsigned> (mxGetN(input[1]))-1;
  vector<PointCloud> scans;
  scans.resize(numScans);
  for (unsigned id = 0; id < numScans; ++id) {
      vector<Surfel3D> *surfels = scans[id].GetPointArray();
      int left = static_cast<int> (offsetData[id]);
      int right = static_cast<int> (offsetData[id+1]);
      surfels->resize(right-left);
      for (int i = left; i < right; ++i) {
          Surfel3D *sur = &(*surfels)[i-left];
          int off = 9*i;
          for (int k = 0; k < 3; ++k) {
              sur->position[k] = pointData[off+k];
              sur->normal[k] = pointData[off+3+k];
              sur->color[k] = 0.f;
          }
      }
  }

  // Setup the parameters 
  double *data = (double*)mxGetData(input[2]);
  float max_distance = static_cast<float> (data[0]);
  float overlap_ratio = static_cast<float> (data[1]);
  float weight_normal = static_cast<float> (data[2]);
  float weight_color = static_cast<float> (data[3]);
  unsigned down_sampling_rate = static_cast<unsigned> (data[4]);
  unsigned num_reweighting_iters = static_cast<unsigned> (data[5]);
  unsigned num_gauss_newton_iters = static_cast<unsigned> (data[6]);
  float weight_point2planeDis = static_cast<float> (data[7]);
  //
  vector<Affine3d> opt_poses;
  opt_poses.resize(numScans);
  vector<PointCorres2> pointcorres;
  
// Perform multi-scan registration
  JointPairwiseRegistration jpr;
  jpr.Compute(scans,
          max_distance, overlap_ratio,
          weight_normal, weight_color,
          down_sampling_rate, num_reweighting_iters, num_gauss_newton_iters,
          weight_point2planeDis,
          &opt_poses,
          &pointcorres);

  // Output optimized scan poses
  output[0] = mxCreateDoubleMatrix(12, numScans, mxREAL);
  data = mxGetPr(output[0]);
  for (unsigned scanid = 0; scanid < numScans; ++scanid) {
	  for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 3; ++j) {
              data[12*scanid+i*3+j] = opt_poses[scanid][i][j];
          }
      }
  }
  
  // Output optimized correspondences
  unsigned numcorres = pointcorres.size();
  output[1] = mxCreateDoubleMatrix(5, numcorres, mxREAL);
  data = mxGetPr(output[1]);
  for (unsigned corid = 0; corid < numcorres; ++corid) {
      const PointCorres2 &pc = pointcorres[corid];
      data[5*corid] = pc.sourceSurfId + 1;
      data[5*corid+1] = pc.targetSurfId + 1;
      data[5*corid+2] = pc.sourcePointId + 1;
      data[5*corid+3] = pc.targetPointId + 1;
      data[5*corid+4] = pc.weight;
  }
}