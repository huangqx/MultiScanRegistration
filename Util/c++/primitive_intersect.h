#ifndef primitive_intersect_h_
#define primitive_intersect_h_

#include "linear_algebra.h"

#define EPSILON 1e-14

struct Ray3D {
 public:
   Ray3D() {}
   Ray3D(const Vector3d &o, const Vector3d &d) {
     origin = o;
     direction = d;
     for (int i = 0; i < 3; ++i) {
       inv_direction[i] = 1/d[i];
       sign[i] = (inv_direction[i] < 0);
     }
   }
   ~Ray3D() {
   }
   Vector3d origin;
   Vector3d direction;
   Vector3d inv_direction;
   int sign[3];
};

struct Box3D {
 public:
   Box3D() {
     corners[0][0] = corners[0][1] = corners[0][2] = 0.0;
     corners[1][0] = corners[1][1] = corners[1][2] = 1.0;
   }
   ~Box3D() {
   }
   Vector3d corners[2];
};

struct IntersectingPoint {
 public:
  IntersectingPoint() {
    faceId = -1;
    paraU = 0.0;
    paraV = 0.0;
    distance = 0.0;
  }
  ~IntersectingPoint() {
  }
  int faceId;
  double paraU;
  double paraV;
  double distance;
};



// Functions that computing the intersection between rays and boxes and triangles
bool TriangleBoxIntersection(const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, const Box3D &box);

bool RayBoxIntersection(const Ray3D &ray, const Box3D &box, const double &maxDis);
bool RayTriangleCullingIntersection(const Ray3D &ray, const Vector3d &v0, const Vector3d &v1, const Vector3d &v2, IntersectingPoint *pt);
bool RayTriangleIntersection(const Ray3D &ray, const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, IntersectingPoint *pt);


#endif