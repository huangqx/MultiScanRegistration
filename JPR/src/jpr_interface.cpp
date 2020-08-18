/*---
function [scanposes, corres] = jpr_interface(datapoints, scan_offsets, paras)
Input: 
Output
--*/

#include "mex.h"
#include "point_cloud.h"
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
  
  

  double *data = (double*)mxGetData(input[2]);
  Vector3d origin;
  origin[0] = data[0];
  origin[1] = data[1];
  origin[2] = data[2];

  double *lookAts = (double*)mxGetData(input[3]);
  unsigned numRays = static_cast<unsigned> (mxGetN(input[3]));

  RayMeshIntersection ray_intersect;
  ray_intersect.Initialize(tri_mesh, 7);

  output[0] = mxCreateDoubleMatrix(4, numRays, mxREAL);
  data = mxGetPr(output[0]);
  
  for (unsigned rayId = 0; rayId < numRays; ++rayId) {
	  Vector3d direction;
	  direction[0] = lookAts[3*rayId];
	  direction[1] = lookAts[3*rayId+1];
	  direction[2] = lookAts[3*rayId+2];
	  direction = direction - origin;
	  direction = direction/sqrt(direction.getSqrNorm());
	  Ray3D ray(origin, direction);

	  IntersectingPoint pt;
	  ray_intersect.Compute(ray, tri_mesh, &pt);
	  data[4*rayId] = pt.faceId + 1;
	  data[4*rayId+1] = pt.paraU;
	  data[4*rayId+2] = pt.paraV;
	  data[4*rayId+3] = pt.distance;
  }
}