#ifndef point_cloud_h_
#define point_cloud_h_

#include "surface.h"

#include <vector>
using	namespace	std;

const unsigned int VERTEX_FLAG_BIT = 1<<8;
const unsigned int VERTEX_ISBOUNDARY_BIT = 1<<9;
const unsigned int VERTEX_ISSEED_BIT = 1<<10;
const unsigned int VERTEX_ISFEATURE_BIT	= 1<<11;
const unsigned int VERTEX_ISACTIVE_BIT = 1<<12;

struct Surfel3D {
public:
  Surfel3D();
  ~Surfel3D();

  int GetMaterial() const;
  void SetMaterial(int material);
  bool GetFlag() const ;
  void SetFlag(bool f);
  bool GetIsFeatureFlag() const;
  void SetIsFeatureFlag(bool f);
  bool GetIsBoundaryFlag() const ;
  void SetIsBoundaryFlag(bool f);
  bool GetIsActiveFlag() const;
  void SetIsActiveFlag(bool f);

  /************************************************************************/
  /* position
  /************************************************************************/
  Vector3f position;

  /************************************************************************/
  /* normal
  /************************************************************************/
  Vector3f normal;

  /************************************************************************/
  /* color 
  /************************************************************************/
  Vector3f color;

  /************************************************************************/
  /* flags
  /************************************************************************/
  int flags;
};

// The following struct stores point-wise correspondences between two point clouds
struct PointCorres2 {
 public:
  PointCorres2() {
    sourceSurfId = targetSurfId = 0;
    sourcePointId = targetPointId = 0;
    weight = 0.f;
  }
  ~PointCorres2() {
  }
  unsigned sourceSurfId;
  unsigned targetSurfId;
  unsigned sourcePointId;
  unsigned targetPointId;
  float weight;
};

class PointCloud: public Surface {
 public:
	PointCloud(){
    current_pose_.T.Initialize();
  }
	virtual ~PointCloud(){
  }
  
  // Bounding box
  void ComputeBoundingBox();

  // Access functions
  vector<Surfel3D>* GetPointArray() {
    return &points_;
  }
  const vector<Surfel3D>* GetPointArray() const {
    return &points_;
  }
 protected:
  vector<Surfel3D> points_; //vertex array
};

#endif