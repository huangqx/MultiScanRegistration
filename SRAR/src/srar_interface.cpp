/*---
function [scanposes, corres] = jpr_interface(datapoints, scan_offsets, paras)
Input: 
Output
--*/

#include "mex.h"
#include "point_cloud.h"
#include "linear_algebra.h"
#include "affine_transformation.h"
#include "simul_reg_and_recons.h"
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
  if (nargout != 3) {
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
              sur->position[k] = static_cast<float> (pointData[off+k]);
              sur->normal[k] = static_cast<float> (pointData[off+3+k]);
              sur->color[k] = 0.f;
          }
      }
  }
  // Compute the bounding box of all the input points
  BoundingBox box_all;
  box_all.Initialize();

  for (unsigned id = 0; id < numScans; ++id) {
    scans[id].ComputeBoundingBox();
    box_all.Insert_A_Box(scans[id].GetBoundingBox());
  }
  double maxDim = max(max(box_all.size[0], box_all.size[1]), box_all.size[2]);

  // Setup the parameters 
  double *data = (double*)mxGetData(input[2]);
  
  // Code for simultaneous registration and reconstruction
  SimulRegAndRecons srar;
  SimulRegAndReconsPara para;
  para.stride = static_cast<int> (data[0]);
  para.weightPoint2PlaneDis = static_cast<double> (data[1]);
  para.minNumPointsPerCell = static_cast<int> (data[2]);
  para.maxNumPointsPerCell = static_cast<int> (data[3]);
  para.numAlternatingIterations = static_cast<int> (data[4]);

  // Perform santity check, so that gridSize_coarse and gridSize_fine fall into
  // suitable ranges (you should set them carefully)
  double gridSize_coarse = static_cast<double> (data[5]);
  double gridSize_fine = static_cast<double> (data[6]);
  gridSize_coarse = max(maxDim / 128.0, min(maxDim / 32.0, gridSize_coarse));
  gridSize_fine = max(maxDim / 384.0, min
          (maxDim / 96.0, gridSize_fine));

  int num_levels = static_cast<int> (data[7]);
  num_levels = min(4, max(1, num_levels));
  
  vector<Affine3d> opt_poses;
  opt_poses.resize(numScans);
  for (int level_id = 0; level_id < num_levels; ++level_id) {
    double t = 0.0;
    if (num_levels > 1)
      t = static_cast<double> (level_id) / (num_levels - 1.0);
    para.gridSize = exp(log(gridSize_coarse) * (1 - t) + log(gridSize_fine) * t);
    //
    srar.AlternatingOpt(scans, para, &opt_poses);
  }
  PointCloud super_scan;
  vector<vector<PointCorres>> corres;
  para.stride = 1;
  srar.LatentSurfOpt(scans, opt_poses, para, &corres, &super_scan);

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
  
  // Output the latent surface
  output[1] = mxCreateDoubleMatrix(9, super_scan.GetPointArray()->size(), mxREAL);
  data = mxGetPr(output[1]);
  for (unsigned pid = 0; pid < super_scan.GetPointArray()->size(); ++pid) {
    Surfel3D *sur = &(*super_scan.GetPointArray())[pid];
    for (int k = 0; k < 3; ++k) {
        data[9*pid+k] = sur->position[k];
        data[9*pid+k+3] = sur->normal[k];
        data[9*pid+k+6] = sur->color[k];
    }
  }
  
  // Output the dense correspondences between the input scans and the latent surface
  unsigned numcorres = 0;
  for (unsigned scanid = 0; scanid < corres.size(); ++scanid) {
      numcorres += corres[scanid].size();
  }
  output[2] = mxCreateDoubleMatrix(4, numcorres, mxREAL);
  data = mxGetPr(output[2]);
  unsigned off = 0;
  for (unsigned scanid = 0; scanid < corres.size(); ++scanid) {
      for (unsigned j = 0; j < corres[scanid].size(); ++j) {
          data[4*off] = scanid + 1;
          data[4*off+1] = corres[scanid][j].sourcePointId + 1;
          data[4*off+2] = corres[scanid][j].targetPointId + 1;
          data[4*off+3] = corres[scanid][j].weight;
          off++;
      }
  }
}