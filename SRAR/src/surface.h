#ifndef surface_h_
#define surface_h_

#include <string>
#include <vector>

#include "bounding_box.h"
#include "linear_algebra.h"
#include "affine_transformation.h"

using namespace std;

const unsigned int SURFACE_RENDER_BIT = 1<<13;
const unsigned int SURFACE_LOCKED_BIT = 1<<14;
const unsigned int SURFACE_ISQUAD_BIT = 1<<15;

struct RigidPose3D {
  RigidPose3D() {
    T.Initialize();
  }
  ~RigidPose3D() {}
  Affine3d T;
  string name;
};

class Surface {
 public:
  Surface();
  virtual ~Surface();
  bool Render() const;
  void SetRenderBit(bool f);
  bool Locked() const;
  void SetLockedBit(bool f);
  bool IsQuad() const;
  void SetIsQuadBit(bool f);
  
  BoundingBox* GetBoundingBox() {
    return &bounding_box_;
  }
  const BoundingBox* GetBoundingBox() const {
    return &bounding_box_;
  }
  RigidPose3D* GetCurrentPose() {
    return &current_pose_;
  }
  const RigidPose3D* GetCurrentPose() const {
    return &current_pose_;
  }

  void SetFileName(wstring &name) {
    filename = name;
  }
  wstring &GetFileName() {
    return filename;
  }
  const wstring &GetFileName() const {
    return filename;
  }
 protected:
  wstring filename;
  BoundingBox	bounding_box_; // Bounding box of the current surface
  vector<RigidPose3D> poses_; // Store all candidate poses used in computation
  RigidPose3D current_pose_; // Current Pose 
 private:
  int flags_; 
};
#endif