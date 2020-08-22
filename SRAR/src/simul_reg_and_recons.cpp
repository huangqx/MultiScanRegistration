#include "simul_reg_and_recons.h"
#include "pairwise_rigid_alignment.h"
#include "bounding_box.h"
#include <algorithm>
#include <set>
using namespace std;

void SimulRegAndRecons::AlternatingOpt(const vector<PointCloud>& scans,
  const SimulRegAndReconsPara& para,
  vector<Affine3d>* opt_poses) {
  if (opt_poses->size() != scans.size())
    return;

  // Allocate space for correspondences between the input scans and the latent surface
  vector<vector<PointCorres>> corres;
  // Allocate space for the latent surface
  PointCloud latentSurf;

  // Parameters used in Pairwise alignment
  PairwiseRigidAlign aligner;
  PairwiseRigidAligPara para_align;
  
  for (int alterOptId = 0; alterOptId < para.numAlternatingIterations; ++alterOptId) {
    // Compute the latent surface and the point-wise correspondences between
    // the input scans and the latent surface
    LatentSurfOpt(scans, *opt_poses,
      para,
      &corres,
      &latentSurf);

    // Perform alignment between each scan and the latent surface
    for (unsigned scanId = 0; scanId < scans.size(); ++scanId) {
      vector<int> movingPointIds;
      vector<int> targetPointIds;
      vector<double> corresWeights;
      const vector<PointCorres>* pt = &corres[scanId];
      movingPointIds.resize(pt->size());
      targetPointIds.resize(pt->size());
      corresWeights.resize(pt->size());
      for (unsigned i = 0; i < pt->size(); ++i) {
        movingPointIds[i] = (*pt)[i].sourcePointId;
        targetPointIds[i] = (*pt)[i].targetPointId;
        corresWeights[i] = (*pt)[i].weight;
      }
      aligner.Compute_KnownCorres(latentSurf, scans[scanId],
        movingPointIds,
        targetPointIds,
        corresWeights,
        para.weightPoint2PlaneDis, 6,
        &(*opt_poses)[scanId]);
    }
  }
}

void SimulRegAndRecons::GenerateSuperScan(const vector<PointCloud>& scans,
  const vector<Affine3d>& opt_poses,
  const SimulRegAndReconsPara& para,
  PointCloud* super_scan) {
  // Create the buffer for all involving points
  unsigned numTotalPoints = 0;
  for (unsigned scanId = 0; scanId < scans.size(); ++scanId)
    numTotalPoints += static_cast<unsigned> (scans[scanId].GetPointArray()->size());

  vector<Surfel3D> all_points;
  all_points.resize(numTotalPoints);
  unsigned off = 0;
  for (unsigned scanId = 0; scanId < scans.size(); ++scanId) {
    const PointCloud& scan = scans[scanId];
    const Affine3d & pose = opt_poses[scanId];
    for (unsigned i = 0; i < scan.GetPointArray()->size(); ++i) {
      const Surfel3D* ori_point = &(*scan.GetPointArray())[i];
      Surfel3D* cur_point = &all_points[off];
      for (int j = 0; j < 3; ++j) {
        cur_point->position[j] = static_cast<float> (pose[0][j]
          + pose[1][j] * ori_point->position[0]
          + pose[2][j] * ori_point->position[1]
          + pose[3][j] * ori_point->position[2]);
        cur_point->normal[j] = static_cast<float> (pose[1][j] * ori_point->normal[0]
          + pose[2][j] * ori_point->normal[1]
          + pose[3][j] * ori_point->normal[2]);
      }
      cur_point->color = ori_point->color;
      off = off + 1;
    }
  }
  Octree3D* octree = GenerateOctree(all_points, para.gridSize, 1);
  if (octree == NULL)
    return;

  // Obtain all leaf nodes
  vector<int> leaf_cell_x_ids, leaf_cell_y_ids, leaf_cell_z_ids;
  octree->CollectAllLeafs(&leaf_cell_x_ids, &leaf_cell_y_ids, &leaf_cell_z_ids);

  int max_cell_coord = (1 << octree->max_depth) - 1;

  // For each leaf cell, perform clustering in its neighbors
  int point_id = 0;
  super_scan->GetPointArray()->resize(leaf_cell_x_ids.size());
  for (unsigned leaf_id = 0; leaf_id < leaf_cell_x_ids.size(); ++leaf_id) {
    const OctreeLeaf3D* center_leaf = octree->QueryLeaf(leaf_cell_x_ids[leaf_id],
      leaf_cell_y_ids[leaf_id],
      leaf_cell_z_ids[leaf_id]);
    if (center_leaf == NULL)
      continue;

    int numpts = static_cast<int> (center_leaf->point_indices.size());
    if (numpts >= para.minNumPointsPerCell) {
      Vector3f pos, nor, color;
      pos[0] = pos[1] = pos[2] = 0.f;
      nor[0] = nor[1] = nor[2] = 0.f;
      color[0] = color[1] = color[2] = 0.f;
      for (unsigned i = 0; i < center_leaf->point_indices.size(); ++i) {
        pos += all_points[center_leaf->point_indices[i]].position;
        nor += all_points[center_leaf->point_indices[i]].normal;
        color += all_points[center_leaf->point_indices[i]].color;
      }
      Surfel3D* pt = &(*super_scan->GetPointArray())[point_id];
      pt->position = pos/static_cast<float> (numpts);
      pt->normal = nor/static_cast<float> (sqrt(nor.getSqrNorm()));
      pt->color = color/static_cast<float> (numpts);
      point_id++;
    }
  }
  super_scan->GetPointArray()->resize(point_id);
  delete octree;
}

void SimulRegAndRecons::SurfelFitting(
  const vector<Vector3f>& poss,
  const vector<Vector3f>& nors,
  const float& weight_nor,
  Surfel3D* fit_sur) {
  unsigned num_points = static_cast<unsigned> (poss.size());
  for (unsigned i = 0; i < 3; ++i) {
    fit_sur->position[i] = 0.f;
    fit_sur->normal[i] = 0.f;
  }
  for (unsigned pid = 0; pid < num_points; ++pid) {
    fit_sur->position += poss[pid];
    fit_sur->normal += nors[pid];
  }
  fit_sur->position /= static_cast<double> (num_points);
  fit_sur->normal.normalize();
}

void SimulRegAndRecons::LatentSurfOpt(const vector<PointCloud>& scans,
  const vector<Affine3d>& opt_poses,
  const SimulRegAndReconsPara& para,
  vector<vector<PointCorres>>* corres,
  PointCloud* latentSurf) {
  unsigned num_total_points = 0;
  for (unsigned id = 0; id < scans.size(); ++id)
    num_total_points += scans[id].GetPointArray()->size();
  
  vector<Surfel3D> all_points;
  vector<unsigned> scan_indices, point_indices;
  all_points.resize(num_total_points);
  scan_indices.resize(num_total_points);
  point_indices.resize(num_total_points);

  unsigned off = 0;
  for (unsigned id = 0; id < scans.size(); ++id) {
    const PointCloud& pc = scans[id];
    const Affine3d& pose = opt_poses[id];
    for (unsigned i = 0; i < pc.GetPointArray()->size(); ++i) {
      const Surfel3D& sur = (*pc.GetPointArray())[i];
      Vector3d cur_pos = pose[0] 
        + pose[1] * sur.position[0] 
        + pose[2] * sur.position[1] 
        + pose[3] * sur.position[2];
      Vector3d cur_nor = pose[1] * sur.normal[0]
        + pose[2] * sur.normal[1]
        + pose[3] * sur.normal[2];
      all_points[off].color = sur.color;
      for (int k = 0; k < 3; ++k) {
        all_points[off].position[k] = cur_pos[k];
        all_points[off].normal[k] = cur_nor[k];
      }
      scan_indices[off] = id;
      point_indices[off] = i;
      off++;
    }
  }

  Octree3D* octree = GenerateOctree(all_points, para.gridSize, para.stride);
  if (octree == NULL)
    return;

  // Obtain all leaf nodes
  vector<int> leaf_cell_x_ids, leaf_cell_y_ids, leaf_cell_z_ids;
  octree->CollectAllLeafs(&leaf_cell_x_ids, &leaf_cell_y_ids, &leaf_cell_z_ids);
  
  // Initialize point correspondences
  vector<unsigned> offsets;
  offsets.resize(scans.size());

  corres->clear();
  corres->resize(scans.size());
  for (unsigned id = 0; id < scans.size(); ++id) {
    offsets[id] = 0;
    (*corres)[id].resize(scans[id].GetPointArray()->size());
  }
  unsigned num_surfels = 0;
  latentSurf->GetPointArray()->resize(leaf_cell_x_ids.size());
  
  for (unsigned leaf_id = 0; leaf_id < leaf_cell_x_ids.size(); ++leaf_id) {
    const OctreeLeaf3D* leaf = octree->QueryLeaf(
      leaf_cell_x_ids[leaf_id],
      leaf_cell_y_ids[leaf_id],
      leaf_cell_z_ids[leaf_id]);

    unsigned num_points = leaf->point_indices.size();
    if (num_points >= para.minNumPointsPerCell) {
      //check whether this cell contains points from two different scans
      unsigned min_scan_id = scans.size(), max_scan_id = 0;
      for (unsigned i = 0; i < leaf->point_indices.size(); ++i) {
        max_scan_id = max(max_scan_id, scan_indices[leaf->point_indices[i]]);
        min_scan_id = min(min_scan_id, scan_indices[leaf->point_indices[i]]);
      }
      if (min_scan_id == max_scan_id)
        continue;

      vector<Vector3f> vertex_poss, vertex_nors;
      vertex_poss.resize(num_points);
      vertex_nors.resize(num_points);
      for (unsigned i = 0; i < leaf->point_indices.size(); ++i) {
        vertex_poss[i] = all_points[leaf->point_indices[i]].position;
        vertex_nors[i] = all_points[leaf->point_indices[i]].normal;
      }
      SurfelFitting(vertex_poss, vertex_nors, 0.0f,
        &(*latentSurf->GetPointArray())[num_surfels]);
      for (unsigned j = 0; j < leaf->point_indices.size(); ++j) {
        unsigned scanId = scan_indices[leaf->point_indices[j]];
        unsigned surfelId = point_indices[leaf->point_indices[j]];
        (*corres)[scanId][offsets[scanId]].sourcePointId = surfelId;
        (*corres)[scanId][offsets[scanId]].targetPointId = num_surfels;
        (*corres)[scanId][offsets[scanId]].weight = 1.f;
        offsets[scanId]++;
      }
      num_surfels++;
    }
  }
  latentSurf->GetPointArray()->resize(num_surfels);
  for (unsigned id = 0; id < corres->size(); ++id)
    (*corres)[id].resize(offsets[id]);
  
  int max_cell_coord = (1 << octree->max_depth) - 1;
  // Release octree
  delete octree;
}

Octree3D* SimulRegAndRecons::GenerateOctree(
  const vector<Surfel3D> &points,
  const double& grid_size,
  const unsigned& stride) {
  // Generate the bounding box
  BoundingBox b_box;
  b_box.Initialize();
  for (unsigned point_id = 0; point_id < points.size(); ++point_id)
    b_box.Insert_A_Point(points[point_id].position);

  double cube_size = max(max(b_box.size[0], b_box.size[1]), b_box.size[2]);
  if (cube_size < grid_size * 2.0)
    return NULL;

  max_depth_ = static_cast<int> (int(log((cube_size / grid_size) + 1.0) / log(2.0))) + 1;
  cube_size = grid_size * (1 << max_depth_);
  Octree3D* octree = new Octree3D();
  octree->max_depth = max_depth_;
  octree->grid_size = static_cast<float> (grid_size);
  octree->left_corner[0] = static_cast<float> (b_box.center_point[0] - (cube_size / 2.0));
  octree->left_corner[1] = static_cast<float> (b_box.center_point[1] - (cube_size / 2.0));
  octree->left_corner[2] = static_cast<float> (b_box.center_point[2] - (cube_size / 2.0));
  octree->root = NULL;

  for (unsigned point_id = 0; point_id < points.size(); ++point_id) {
    const Vector3f& pos = points[point_id].position;
    cell_offsets_[0] = static_cast<int> ((pos[0] - octree->left_corner[0]) / grid_size);
    cell_offsets_[1] = static_cast<int> ((pos[1] - octree->left_corner[1]) / grid_size);
    cell_offsets_[2] = static_cast<int> ((pos[2] - octree->left_corner[2]) / grid_size);
    if (cell_offsets_[0] % stride != 0
      || cell_offsets_[1] % stride != 0
      || cell_offsets_[2] % stride != 0)
      continue;

    depth_ = 0;
    InsertAPoint(point_id, &octree->root);
  }
  return octree;
}

void SimulRegAndRecons::InsertAPoint(const int& point_index,
  Node3D** current_node) {
  if (depth_ < max_depth_) {
    if (*current_node == NULL)
      * current_node = new OctreeNode3D();

    OctreeNode3D* nl_node = (OctreeNode3D*)(*current_node);
    int half_box_width = 1 << (max_depth_ - depth_ - 1);
    int i = 0, j = 0, k = 0;
    if (cell_offsets_[0] >= half_box_width) {
      cell_offsets_[0] -= half_box_width;
      i = 1;
    }
    if (cell_offsets_[1] >= half_box_width) {
      cell_offsets_[1] -= half_box_width;
      j = 1;
    }
    if (cell_offsets_[2] >= half_box_width) {
      cell_offsets_[2] -= half_box_width;
      k = 1;
    }
    int child_offset = (i << 2) + (j << 1) + k;
    depth_++;
    InsertAPoint(point_index, &(nl_node->children[child_offset]));
  }
  else {
    if (*current_node == NULL)
      * current_node = new OctreeLeaf3D();

    OctreeLeaf3D* ol_node = (OctreeLeaf3D*)(*current_node);
    ol_node->point_indices.push_back(point_index);
  }
}