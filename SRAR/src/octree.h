#ifndef octree_h_
#define octree_h_

#include "linear_algebra.h"

#include <vector>
using namespace std;

struct Node3D {
public:
  virtual ~Node3D() {
  }
};

// Regular octree node 
struct OctreeNode3D : public Node3D {
public:
  OctreeNode3D() {
    for (int i = 0; i < 8; i++) {
      children[i] = NULL;
    }
  }
  ~OctreeNode3D() {
    for (int i = 0; i < 8; i++) {
      if (children[i] != NULL) {
        delete children[i];
        children[i] = NULL;
      }
    }
  }
  Node3D* children[8];
};

// octree leaf node
struct OctreeLeaf3D : public Node3D {
public:
  OctreeLeaf3D() {
  }
  ~OctreeLeaf3D() {
  }

  vector<int> point_indices;
};

// octree struct which is also regarded as a compact representation of a sub-grid
struct Octree3D {
public:
  Octree3D() {
    root = NULL;
  }
  ~Octree3D() {
    if (root != NULL) {
      delete root;
      root = NULL;
    }
  }
  // Query operations
  const OctreeLeaf3D *QueryLeaf(int cell_x_id, int cell_y_id,
    int cell_z_id) const;
  void CollectAllLeafs(vector<int> *cell_x_ids, vector<int> *cell_y_ids,
    vector<int> *cell_z_ids) const;

  Node3D *root;

  Vector3f left_corner;
  int max_depth;
  float grid_size;
private:
  const OctreeLeaf3D *Iterate1(Node3D *cur_node) const ;
  void Iterate2(Node3D *cur_node,
    vector<int> *cell_x_ids,
    vector<int> *cell_y_ids,
    vector<int> *cell_z_ids ) const;
}; 

#endif