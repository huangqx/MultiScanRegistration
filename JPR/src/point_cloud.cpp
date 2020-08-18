#include "point_cloud.h"

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

void PointCloud::Read(FILE *file_ptr) {
  ReadPoints(file_ptr);
  bounding_box_.Read(file_ptr);
  fread(&current_pose_.T, sizeof(Affine3d), 1, file_ptr);
}

void PointCloud::Write(FILE *file_ptr) {
  WritePoints(file_ptr);
  bounding_box_.Write(file_ptr);
  fwrite(&current_pose_.T, sizeof(Affine3d), 1, file_ptr);
}

void PointCloud::ReadPoints(FILE *file_ptr) {
  int	length;
  fread(&length, sizeof(int), 1, file_ptr);
  points_.resize(length);
  if (length > 0)
    fread(&points_[0], sizeof(Surfel3D)*length, 1, file_ptr);
}

void	PointCloud::WritePoints(FILE *file_ptr) {
  int	length = static_cast<int> (points_.size());
  fwrite(&length, sizeof(int), 1, file_ptr);
  if (length > 0)
    fwrite(&points_[0], sizeof(Surfel3D)*length, 1, file_ptr);
}

void SuperPointCloud::Read(FILE* file_ptr) {
  pc_.Read(file_ptr);
  int	length;
  fread(&length, sizeof(int), 1, file_ptr);
  child_pc_ids_.resize(length);
  if (length > 0)
    fread(&child_pc_ids_[0], sizeof(unsigned) * length, 1, file_ptr);
  fread(&length, sizeof(int), 1, file_ptr);
  child_pc_poses_.resize(length);
  if (length > 0)
    fread(&child_pc_poses_[0], sizeof(Affine3d) * length, 1, file_ptr);
}

void SuperPointCloud::Write(FILE* file_ptr) {
  pc_.Write(file_ptr);
  int	length = static_cast<int> (child_pc_ids_.size());
  fwrite(&length, sizeof(int), 1, file_ptr);
  if (length > 0)
    fwrite(&child_pc_ids_[0], sizeof(unsigned) * length, 1, file_ptr);
  length = static_cast<int> (child_pc_poses_.size());
  fwrite(&length, sizeof(int), 1, file_ptr);
  if (length > 0)
    fwrite(&child_pc_poses_[0], sizeof(Affine3d) * length, 1, file_ptr);
}