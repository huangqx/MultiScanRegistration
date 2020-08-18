#include "octree.h"
#include <algorithm>
#include <set>
using namespace std;

int depth_;
int cell_offsets_[3];

const OctreeLeaf3D *Octree3D::QueryLeaf(int cell_x_id, int cell_y_id,
  int cell_z_id) const {
  depth_ = 0;
  cell_offsets_[0] = cell_x_id;
  cell_offsets_[1] = cell_y_id;
  cell_offsets_[2] = cell_z_id;
  return Iterate1(root);
}
void Octree3D::CollectAllLeafs(vector<int> *cell_x_ids,
  vector<int> *cell_y_ids,
  vector<int> *cell_z_ids) const {
  depth_ = 0;
  cell_x_ids->clear();
  cell_y_ids->clear();
  cell_z_ids->clear();
  cell_offsets_[0] = cell_offsets_[1] = cell_offsets_[2] = 0;
  Iterate2(root, cell_x_ids, cell_y_ids, cell_z_ids);
}

const OctreeLeaf3D* Octree3D::Iterate1(Node3D *cur_node) const {
  if (depth_ < max_depth) {
    OctreeNode3D *nl_node = (OctreeNode3D*)cur_node;
    int half_box_width = 1<<(max_depth - depth_ - 1);
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
    int child_offset = (i<<2) + (j<<1) + k;
    if (nl_node->children[child_offset] == NULL)
      return NULL;
    else {
      depth_++;
      return Iterate1(nl_node->children[child_offset]);
    }
  } else {
    if (cur_node == NULL)
      return NULL;
    else
      return (OctreeLeaf3D*)cur_node;
  }
  return NULL;
}

void Octree3D::Iterate2(Node3D *cur_node,
  vector<int> *cell_x_ids, vector<int> *cell_y_ids, vector<int> *cell_z_ids) const {
  if (depth_ < max_depth) {
    OctreeNode3D *nl_node = (OctreeNode3D*)cur_node;
    int current_depth = depth_;
    int old_cell_offset[3];
    old_cell_offset[0] = cell_offsets_[0];
    old_cell_offset[1] = cell_offsets_[1];
    old_cell_offset[2] = cell_offsets_[2];
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          int child_offset = (i<<2) + (j<<1) + k;
          if (nl_node->children[child_offset] != NULL) {
            cell_offsets_[0] = old_cell_offset[0] + (i<<(max_depth-current_depth-1));
            cell_offsets_[1] = old_cell_offset[1] + (j<<(max_depth-current_depth-1));
            cell_offsets_[2] = old_cell_offset[2] + (k<<(max_depth-current_depth-1));
            depth_ = current_depth + 1;
            Iterate2(nl_node->children[child_offset],
              cell_x_ids, cell_y_ids, cell_z_ids);
          }
        }
      }
    }
  } else {
    cell_x_ids->push_back(cell_offsets_[0]);
    cell_y_ids->push_back(cell_offsets_[1]);
    cell_z_ids->push_back(cell_offsets_[2]);
  }
}

