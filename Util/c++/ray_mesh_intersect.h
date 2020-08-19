#ifndef ray_mesh_intersect_h_
#define ray_mesh_intersect_h_

#include "primitive_intersect.h"
#include "tri_mesh.h"
#include <string>
#include <vector>
using namespace std;


// Virtual base node
class Node3D {
 public:
	virtual ~Node3D() {
	}
//	virtual void queryNode(const Ray3D &Ray, const TriMesh &mesh, IntersectingPoint *pt) = 0;
};

// Non-leaf node in the octree
class OctreeNode3D : public Node3D {
 public:
  OctreeNode3D() {
    for (int i = 0; i < 8; ++i)
      m_children[i] = NULL;
  }
	virtual ~OctreeNode3D()  {
    for (int i = 0; i < 8; ++i) {
      if (m_children[i] != NULL) {
        delete m_children[i];
        m_children[i] = NULL;
      }
    }
	}
	Node3D*			m_children[8];		
	
//  void queryNode(const Ray3D &Ray, const TriMesh &mesh, IntersectingPoint *pt);
};

// Leaf node in the octree
class OctreeLeaf3D : public Node3D {
public:
	/**
	 * the primitives of this leaf
	 */
  vector<unsigned> intersectingFaceIndices;

//	void queryNode(const Ray3D &Ray, const TriMesh &mesh, IntersectingPoint *pt);
};

struct MeshPoint {
 public:
  MeshPoint() {
    pixelId = 0;
    faceId = 0;
    paraU = 0.f;
    paraV = 0.f;
  }
  ~MeshPoint() {
  }
  int pixelId;
  int faceId;
  float paraU;
  float paraV;
};

struct Camera3D {
 public:
  Camera3D() {
    Initialize_1920_1080();
  }
  ~Camera3D() {
  }
  void Initialize_1920_1080() {
    origin[0] = origin[1] = origin[2] = 0;
    origin[2] = 4.0;
    viewDistance = 4.0;

    frontVec[0] = frontVec[1] = frontVec[2] = 0;
    frontVec[2] = -1.0;

    upVec[0] = upVec[1] = upVec[2] = 0;
    upVec[1] = 1.0;

    rightVec[0] = rightVec[1] = rightVec[2] = 0;
    rightVec[0] = 1.0;

    dimRow = 1080;
    dimCol = 1920;
    pixelOffCol = 960;
    pixelOffRow = 540;
    scale = 1080;
  }
  // coordinate system
  Vector3d origin;
  Vector3d frontVec;
  Vector3d upVec;
  Vector3d rightVec;
  double viewDistance;

  double pixelOffRow;
  double pixelOffCol;
  double scale;

  unsigned dimRow;
  unsigned dimCol;

  string name;
};

class RayMeshIntersection {
 public:
   RayMeshIntersection() {
     root_ = NULL;
     numVisitedFaces_ = 0;
   }
   ~RayMeshIntersection() {
     if (root_ != NULL) {
       delete root_;
       root_ = NULL;
     }
   }
   // If intersecting, return the intersecting point
   // O + tD = (1-u-v)*V0 + u*V1 + v*V2
   void Initialize(const TriMesh &mesh, const unsigned &maxDepth);
   void CameraSimulation(const Camera3D &camera, const TriMesh &mesh, vector<MeshPoint> *intersectingPixels);
   bool Compute(const Ray3D &ray, const TriMesh &mesh, IntersectingPoint *pt);
   bool Intersect_BF(const Ray3D &ray, const TriMesh &mesh, IntersectingPoint *pt);
private:
  void CreateTree(const unsigned &faceId,
     const Vector3f &v0, const Vector3f &v1, const Vector3f &v2,
     OctreeNode3D *node);
  void queryNode(const Node3D *curNode, const Ray3D &Ray, const TriMesh &mesh, IntersectingPoint *pt);
private:
  // root of the tree
  OctreeNode3D *root_;

  // maximum depth of the octree
  unsigned maxDepth_;

  // variables that capture the bounding box of the octree
  Vector3d lowerCorner_;
  Vector3d leafNodeSize_;

  // Store the indices of traversed faces
  unsigned numVisitedFaces_;
  vector<int> visitedFaceIndices_;
  vector<bool> visitedFaceFlags_;


  unsigned lowerCornerOff_X_;
  unsigned lowerCornerOff_Y_;
  unsigned lowerCornerOff_Z_;
  int dis2Leaf_;
};

#endif