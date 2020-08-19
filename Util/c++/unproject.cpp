/*---
function [intersectingPoints] = unproject(mesh.vertexPoss, mesh.faceVIds, origin, lookAts)
Input: 

Output

--*/

#include "mex.h"
#include "tri_mesh.h"
#include "ray_mesh_intersect.h"
#include "primitive_intersect.h"
#include "linear_algebra.h"
#include "math.h"
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
  if (nargin != 4) {
    mexErrMsgTxt("Three input arguments required.");
  }
  if (nargout != 1) {
    mexErrMsgTxt("Incorrect number of output arguments."); 
  }

  double *vertexData = (double*)mxGetData(input[0]);
  unsigned numVertices = static_cast<unsigned> (mxGetN(input[0]));
  double *faceData = (double*)mxGetData(input[1]);
  unsigned numFaces = static_cast<unsigned> (mxGetN(input[1]));

  TriMesh tri_mesh;
  tri_mesh.InitializeFromMatlab(vertexData, numVertices,
	  faceData, numFaces);

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


