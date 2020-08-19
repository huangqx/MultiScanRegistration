#include "primitive_intersect.h"
#include "math.h"

bool RayBoxIntersection(const Ray3D &ray, const Box3D &box, const double &maxDis) {
  double tmin = (box.corners[ray.sign[0]][0] - ray.origin[0]) * ray.inv_direction[0];
  double tmax = (box.corners[1-ray.sign[0]][0] - ray.origin[0]) * ray.inv_direction[0];
  double tymin = (box.corners[ray.sign[1]][1] - ray.origin[1]) * ray.inv_direction[1];
  double tymax = (box.corners[1-ray.sign[1]][1] - ray.origin[1]) * ray.inv_direction[1];

  if ( (tmin > tymax) || (tymin > tmax) )
    return false;
  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;

  double tzmin = (box.corners[ray.sign[2]][2] - ray.origin[2]) * ray.inv_direction[2];
  double tzmax = (box.corners[1-ray.sign[2]][2] - ray.origin[2]) * ray.inv_direction[2];

  if ( (tmin > tzmax) || (tzmin > tmax) )
    return false;
  if (tzmin > tmin)
    tmin = tzmin;
  if (tzmax < tmax)
    tmax = tzmax;

  return tmin < maxDis;
  /*  } else {
  Vector3f verts[8];
  for (int i = 0; i < 2; ++i) {
  for (int j = 0; j < 2; ++j) {
  for (int k = 0; k < 2; ++k) {
  int nodeId = 4*i + 2*j + k;
  verts[nodeId][0] = box.corners[i][0];
  verts[nodeId][1] = box.corners[j][1];
  verts[nodeId][2] = box.corners[k][2];
  }
  }
  }
  IntersectingPoint pt;
  pt.distance  = 1e10;
  if (RayTriangleIntersection(ray, verts[0], verts[1], verts[3], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[0], verts[2], verts[3], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[4], verts[5], verts[7], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[4], verts[6], verts[7], &pt))
  return true;

  if (RayTriangleIntersection(ray, verts[0], verts[2], verts[6], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[0], verts[4], verts[6], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[1], verts[3], verts[7], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[1], verts[5], verts[7], &pt))
  return true;

  if (RayTriangleIntersection(ray, verts[0], verts[1], verts[5], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[0], verts[4], verts[5], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[2], verts[3], verts[7], &pt))
  return true;
  if (RayTriangleIntersection(ray, verts[2], verts[6], verts[7], &pt))
  return true;

  return false;
  }*/
}

bool RayTriangleCullingIntersection(const Ray3D &ray,
  const Vector3d&v0, const Vector3d &v1, const Vector3d &v2,
  IntersectingPoint *pt) {
    Vector3d edge1 = v1 - v0;
    Vector3d edge2 = v2 - v0;

    Vector3d pvec = ray.direction.crossProduct(edge2);

    double det = edge1*pvec;

    if (det < EPSILON)
      return false;

    Vector3d tvec = ray.origin - v0;

    // Calculate U parameter
    pt->paraU = tvec*pvec;
    if (pt->paraU < 0.f || pt->paraU > det)
      return false;

    // Calculate V parameter
    Vector3d qvec = tvec.crossProduct(edge1);
    pt->paraV = ray.direction*qvec;
    if (pt->paraV < 0.f || pt->paraU + pt->paraV > det)
      return 0.f;

    // Calculate distance
    pt->distance = edge2*qvec;

    double inv_det = 1.0/det;
    pt->distance *= inv_det;
    pt->paraU *= inv_det;
    pt->paraV *= inv_det;

    return true;
}

bool RayTriangleIntersection(const Ray3D &ray,
  const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, IntersectingPoint *pt) {
  Vector3d edge1, edge2, tvec;
    
  for (int i = 0; i < 3; ++i) {
    edge1[i] = v1[i] - v0[i];
    edge2[i] = v2[i] - v0[i];
    tvec[i] = ray.origin[i] - v0[i];
 }
  Vector3d pvec = ray.direction.crossProduct(edge2);

  double det = edge1*pvec;

  if (det > -EPSILON && det < EPSILON)
    return false;

  double inv_det = 1.0/det;

  /* prepare to test V parameter*/
  Vector3d qvec = tvec.crossProduct(edge1);
  double distance = (edge2*qvec)*inv_det;
  if (distance >= pt->distance)
    return false;

  // calculate U parameter and test bounds
  double paraU = (tvec*pvec)*inv_det;
  if (paraU < 0.0 - EPSILON || paraU > 1.0 + EPSILON)
    return false;

  // Calculate V parameter and test bounds
  double paraV = (ray.direction*qvec)*inv_det;
  if (paraV < 0.0 - EPSILON || paraU + paraV > 1.0 + EPSILON)
    return false;

  pt->distance = distance;
  pt->paraU = paraU;
  pt->paraV = paraV;
  return true;
}


/************************************************************************/
/* The following routines are used in triangle box intersection
/************************************************************************/
bool AXISTEST_X01(const double &a, const double &b, const double &fa, const double &fb,
  const Vector3d &v0, const Vector3d &v2, const Vector3d &maxbox) {
  double minVal = a*v0[1] - b*v0[2];	
  double maxVal= a*v2[1] - b*v2[2];	
  if (minVal > maxVal) {
    double tmp=maxVal;
    maxVal=minVal;
    minVal=tmp;	
  }	
  double rad = fa * maxbox[1] + fb * maxbox[2];
  if (minVal > rad || maxVal < -rad)
    return false;

  return true;
}

  
//! TO BE DOCUMENTED
bool AXISTEST_X2(const double &a, const double &b, const double &fa, const double &fb,
  const Vector3d &v0, const Vector3d &v1, const Vector3d &maxbox)	{
  double minVal = a*v0[1] - b*v0[2];	
  double maxVal = a*v1[1] - b*v1[2];
  if (minVal > maxVal) {
    double tmp = maxVal;
    maxVal = minVal;
    minVal = tmp;	
  }	
  double rad = fa * maxbox[1] + fb * maxbox[2];
  if (minVal > rad || maxVal < -rad)
    return false;
  return true;
}

//! TO BE DOCUMENTED
bool AXISTEST_Y02(const double &a, const double &b, const double &fa, const double &fb,
  const Vector3d &v0, const Vector3d &v2, const Vector3d &maxbox)	{
  double minVal = b*v0[2] - a*v0[0];
  double maxVal = b*v2[2] - a*v2[0];
  if (minVal > maxVal) {
    double tmp = maxVal;
    maxVal = minVal;
    minVal = tmp;
  }
  double rad = fa * maxbox[0] + fb * maxbox[2];
  if (minVal > rad || maxVal < -rad)
    return false;
  return true;
}

//! TO BE DOCUMENTED
bool AXISTEST_Y1(const double &a, const double &b, const double &fa, const double &fb,
  const Vector3d &v0, const Vector3d &v1, const Vector3d &maxbox) {
  double minVal = b*v0[2] - a*v0[0];
  double maxVal = b*v1[2] - a*v1[0];
  
  if (minVal>maxVal) {
    double tmp = maxVal;
    maxVal = minVal;
    minVal = tmp;
  }	
  double rad = fa * maxbox[0] + fb * maxbox[2];	
  if (minVal>rad || maxVal<-rad)
    return false;
  
  return true;
}

//! TO BE DOCUMENTED
bool AXISTEST_Z12(const double &a, const double &b, const double &fa, const double &fb,
  const Vector3d &v1, const Vector3d &v2, const Vector3d &maxbox)	 {
  double minVal = a*v1[0] - b*v1[1];
  double maxVal = a*v2[0] - b*v2[1];
  if (minVal>maxVal) {
    double tmp = maxVal;
    maxVal = minVal;
    minVal = tmp;	
  }	
  
  double rad = fa * maxbox[0] + fb * maxbox[1];
  if (minVal>rad || maxVal<-rad)
    return false;
 
  return true;
}

//! TO BE DOCUMENTED
bool AXISTEST_Z0(const double &a, const double &b, const double &fa, const double &fb,
  const Vector3d &v0, const Vector3d &v1, const Vector3d &maxbox) {
  double minVal = a*v0[0] - b*v0[1];
  double maxVal = a*v1[0] - b*v1[1];
  if (minVal>maxVal) {
    double tmp = maxVal;
    maxVal = minVal;
    minVal = tmp;	
  }	
  double rad = fa * maxbox[0] + fb * maxbox[1];	
  if (minVal>rad || maxVal<-rad) 
    return false;
  
  return true;
}

bool planeBoxOverlap(const Vector3d& normal, const double &dis, const Vector3d& maxbox) {
  Vector3d vmin, vmax;
  for (int q=0;q<=2;q++) {
    if (normal[q] > 0.0)	{
      vmin[q]=-maxbox[q];
      vmax[q]=maxbox[q]; 
    } else {
      vmin[q]=maxbox[q];
      vmax[q]=-maxbox[q];
    }
  }
  if ((normal*vmin)+dis > 0.0) 
    return false;
  if ((normal*vmax)+dis >=0.0) 
    return true;

  return false;
}

bool InsideBox(const Box3D &box, const Vector3f &point) {
  for (int i = 0; i < 3; ++i) {
    if (point[i] < box.corners[0][i] || point[i] > box.corners[1][i])
      return false;
  }
  return true;
}

bool TriangleBoxIntersection(const Vector3f &tv0, const Vector3f &tv1, const Vector3f &tv2, const Box3D &box1) {
  if (InsideBox(box1, tv0) || InsideBox(box1, tv1) || InsideBox(box1, tv2))
    return true;

  Vector3d center = (box1.corners[0] + box1.corners[1])/2;
  Vector3d maxbox = (box1.corners[1] - box1.corners[0])/2;

  Vector3d v0, v1, v2;
  for (int i = 0; i < 3; ++i) {
    v0[i] = tv0[i] - center[i];
    v1[i] = tv1[i] - center[i];
    v2[i] = tv2[i] - center[i];
  }
  Vector3d e0 = v1 - v0;
  Vector3d e1 = v2 - v1;
  Vector3d e2 = v0 - v2;
  Vector3d normal = e0.crossProduct(e1);
  double dis = -normal*v0;

  if (!planeBoxOverlap(normal, dis, maxbox))
    return false;

  // type 3 tests, there are 9 of them

  double fex = fabs(e0[0]);
  double fey = fabs(e0[1]);
  double fez = fabs(e0[2]);
  if (!AXISTEST_X01(e0[2], e0[1], fez, fey, v0, v2, maxbox))
    return false;
  if (!AXISTEST_Y02(e0[2], e0[0], fez, fex, v0, v2, maxbox))
    return false;
  if (!AXISTEST_Z12(e0[1], e0[0], fey, fex, v1, v2, maxbox))
    return false;

  fex = fabs(e1[0]);
  fey = fabs(e1[1]);
  fez = fabs(e1[2]);

  if (!AXISTEST_X01(e1[2], e1[1], fez, fey, v0, v2, maxbox))
    return false;
  if (!AXISTEST_Y02(e1[2], e1[0], fez, fex, v0, v2, maxbox))
    return false;
  if (!AXISTEST_Z0(e1[1], e1[0], fey, fex, v0, v1, maxbox))
    return false;

  fex = fabs(e2[0]);
  fey = fabs(e2[1]);
  fez = fabs(e2[2]);
  if (!AXISTEST_X2(e2[2], e2[1], fez, fey, v0, v1, maxbox))
    return false;
  if (!AXISTEST_Y1(e2[2], e2[0], fez, fex, v0, v1, maxbox))
    return false;
  if (!AXISTEST_Z12(e2[1], e2[0], fey, fex, v1, v2, maxbox))
    return false;

  return true;
}

