#include "point_cloud.h"
#include "linear_algebra.h"

Surfel3D::Surfel3D() {
  flags = 0;
  position[0] = position[1] = position[2] = 0.f;
  normal[0] = normal[1] = normal[2] = 0.f;
  color[0] = color[1] = color[2] = 0.f;
}

Surfel3D::~Surfel3D() {
}

int Surfel3D::GetMaterial() const{
  return flags & 255;
}

void Surfel3D::SetMaterial(int material) {
  material = material & 255;
  flags = flags | 255;
  flags -= 255 - material;
}

bool Surfel3D::GetFlag() const{
  return (flags & VERTEX_FLAG_BIT) == VERTEX_FLAG_BIT;
}

void Surfel3D::SetFlag(bool f) {
  if(f) {
    flags = flags | VERTEX_FLAG_BIT;
  } else {
    flags = flags & (~VERTEX_FLAG_BIT);
  }
}

bool Surfel3D::GetIsBoundaryFlag() const{
  return (flags & VERTEX_ISBOUNDARY_BIT) == VERTEX_ISBOUNDARY_BIT;
}

void Surfel3D::SetIsBoundaryFlag(bool f) {
  if(f) {
    flags = flags | VERTEX_ISBOUNDARY_BIT;
  } else {
    flags = flags & (~VERTEX_ISBOUNDARY_BIT);
  }
} 

bool Surfel3D::GetIsFeatureFlag() const {
  return (flags & VERTEX_ISFEATURE_BIT) == VERTEX_ISFEATURE_BIT;
}

void Surfel3D::SetIsFeatureFlag(bool f) {
  if (f) {
    flags = flags | VERTEX_ISFEATURE_BIT;
  } else {
    flags = flags & (~VERTEX_ISFEATURE_BIT);
  }
}

bool Surfel3D::GetIsActiveFlag() const {
  return (flags & VERTEX_ISACTIVE_BIT) == VERTEX_ISACTIVE_BIT;
}

void Surfel3D::SetIsActiveFlag(bool f) {
  if (f) {
    flags = flags | VERTEX_ISACTIVE_BIT;
  } else {
    flags = flags & (~VERTEX_ISACTIVE_BIT);
  }
}

void PointCloud::ComputeBoundingBox() {
  bounding_box_.Initialize();
  for (unsigned v_id = 0; v_id < points_.size(); ++v_id) {
    Surfel3D &point = points_[v_id];
    Vector3d pos = current_pose_.T[0]
      + current_pose_.T[1] * point.position[0]
      + current_pose_.T[2] * point.position[1]
      + current_pose_.T[3] * point.position[2];
    bounding_box_.Insert_A_Point(point.position);
  }
}