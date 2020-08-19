#include "ray_mesh_intersect.h"
#include "math.h"

// These variables are set and fixed
//Vector3d lowerCornerPos;
//Vector3d  leafNodeSize;

void RayMeshIntersection::CameraSimulation(const Camera3D &camera, const TriMesh &mesh, vector<MeshPoint> *intersectingPixels) {
  intersectingPixels->clear();
 
  unsigned numAlloc = 16384, numIntersects = 0;
  intersectingPixels->resize(numAlloc);


  for (unsigned colId = 0; colId < camera.dimCol; ++colId) {
    for (unsigned rowId = 0; rowId < camera.dimRow; ++rowId) {
      int pixelId = static_cast<int> (colId*camera.dimRow + rowId);
      // Generate the ray in the camera coordinate system
      double coordX = (colId + 0.5 - camera.pixelOffCol)/camera.scale;
      double coordY = -(rowId + 0.5 - camera.pixelOffRow)/camera.scale;
      double coordZ = 1.0;

      Vector3d direc = camera.rightVec*coordX + camera.upVec*coordY + camera.frontVec*coordZ;
      direc.normalize();
      Ray3D ray(camera.origin, direc);
      IntersectingPoint new_pt;

      if (Compute(ray, mesh, &new_pt)) {
        if (numIntersects + 1 > numAlloc) {
          intersectingPixels->resize(2*numAlloc);
          numAlloc *= 2;
        }
        MeshPoint *mp = &(*intersectingPixels)[numIntersects];
        mp->pixelId = pixelId;
        mp->faceId = new_pt.faceId;
        mp->paraU = static_cast<float> (new_pt.paraU);
        mp->paraV = static_cast<float> (new_pt.paraV);

        numIntersects++;
      }
    }
  }
  intersectingPixels->resize(numIntersects);
}

bool RayMeshIntersection::Compute(const Ray3D &ray, const TriMesh &mesh, IntersectingPoint *pt) {
  lowerCornerOff_X_ = 0;
  lowerCornerOff_Y_ = 0;
  lowerCornerOff_Z_ = 0;

  dis2Leaf_ = maxDepth_;
//  lowerCornerPos = lowerCorner_;
//  leafNodeSize = leafNodeSize_;

  pt->faceId = -1;
  pt->distance = 1e10;
//  root_->queryNode(ray, mesh, pt);
  queryNode(root_, ray, mesh, pt);

  for (unsigned tryId = 0; tryId < numVisitedFaces_; ++tryId)
    visitedFaceFlags_[visitedFaceIndices_[tryId]] = false;

  numVisitedFaces_ = 0;

  if (pt->faceId == -1)
    return false;
  else
    return true;
}

void RayMeshIntersection::queryNode(const Node3D *curNode, const Ray3D &ray, const TriMesh &mesh, IntersectingPoint *pt) {
  int width = 1<<dis2Leaf_;

  Box3D box;
  box.corners[0][0] = lowerCorner_[0] + lowerCornerOff_X_*leafNodeSize_[0];
  box.corners[0][1] = lowerCorner_[1] + lowerCornerOff_Y_*leafNodeSize_[1];
  box.corners[0][2] = lowerCorner_[2] + lowerCornerOff_Z_*leafNodeSize_[2];
  box.corners[1][0] = lowerCorner_[0] + (lowerCornerOff_X_+width)*leafNodeSize_[0];
  box.corners[1][1] = lowerCorner_[1] + (lowerCornerOff_Y_+width)*leafNodeSize_[1];
  box.corners[1][2] = lowerCorner_[2] + (lowerCornerOff_Z_+width)*leafNodeSize_[2];

  if (!RayBoxIntersection(ray, box, pt->distance)) {
    return;
  }

  if (dis2Leaf_ > 0) {
    unsigned old_lowerCornerOff_X = lowerCornerOff_X_;
    unsigned old_lowerCornerOff_Y = lowerCornerOff_Y_;
    unsigned old_lowerCornerOff_Z = lowerCornerOff_Z_;
    int old_dis2Leaf = dis2Leaf_;

    const OctreeNode3D *node = (OctreeNode3D*)curNode;

    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k  = 0; k < 2; ++k) {
          unsigned nodeId = 4*i + 2*j + k;
          if (node->m_children[nodeId] == NULL)
            continue;
          else {
            dis2Leaf_ = old_dis2Leaf - 1;
            lowerCornerOff_X_ = old_lowerCornerOff_X + i*(width/2);
            lowerCornerOff_Y_ = old_lowerCornerOff_Y + j*(width/2);
            lowerCornerOff_Z_ = old_lowerCornerOff_Z + k*(width/2);
            queryNode(node->m_children[nodeId], ray, mesh, pt);
          }
        }
      }
    }
  } else {
    const OctreeLeaf3D *leaf = (OctreeLeaf3D*)curNode;

    // check each triangle separately
    // intersectingFaceIndices
    for (unsigned id = 0; id < leaf->intersectingFaceIndices.size(); ++id) {
      unsigned faceId = leaf->intersectingFaceIndices[id];
      if (visitedFaceFlags_[faceId])
        continue;
      const TriFace &face = mesh.faces[faceId];
      const Vector3f &v0 = mesh.vertices[face.indices[0]].pos;
      const Vector3f &v1 = mesh.vertices[face.indices[1]].pos;
      const Vector3f &v2 = mesh.vertices[face.indices[2]].pos;

      if (RayTriangleIntersection(ray, v0, v1, v2, pt))
        pt->faceId = faceId;

      visitedFaceIndices_[numVisitedFaces_] = faceId;
      visitedFaceFlags_[faceId] = true;
      numVisitedFaces_++;
    }
  }
}

bool RayMeshIntersection::Intersect_BF(const Ray3D &ray, const TriMesh &mesh, IntersectingPoint *pt) {
  pt->faceId = -1;
  pt->distance = 1e10;
  for (unsigned fId = 0; fId < mesh.faces.size(); ++fId) {
    const TriFace &face = mesh.faces[fId];
    const Vector3f &v0 = mesh.vertices[face.indices[0]].pos;
    const Vector3f &v1 = mesh.vertices[face.indices[1]].pos;
    const Vector3f &v2 = mesh.vertices[face.indices[2]].pos;


    if (RayTriangleIntersection(ray, v0, v1, v2, pt))
      pt->faceId = fId;
  }
  return pt->faceId >= 0;
}
void RayMeshIntersection::Initialize(const TriMesh &mesh, const unsigned &maxDepth) {
  maxDepth_ = maxDepth;
  BoundingBox bbox;
  bbox.Initialize();
  for (unsigned vId = 0; vId < mesh.vertices.size(); ++vId)
    bbox.Insert_A_Point(mesh.vertices[vId].pos);

  double d = bbox.size[0];
  if (d < bbox.size[1])
      d = bbox.size[1];
  if (d < bbox.size[2])
      d = bbox.size[2];
  d = d + 1e-4;
  
  leafNodeSize_[0] = d/(1<<maxDepth_);
  leafNodeSize_[1] = d/(1<<maxDepth_);
  leafNodeSize_[2] = d/(1<<maxDepth_);

  lowerCorner_[0] = bbox.center_point[0] - d/2;
  lowerCorner_[1] = bbox.center_point[1] - d/2;
  lowerCorner_[2] = bbox.center_point[2] - d/2;

  if (root_ != NULL) {
    delete root_;
    root_ = NULL;
  }

  root_ = new OctreeNode3D();

//  lowerCornerPos = lowerCorner_;
//  leafNodeSize = leafNodeSize_;
  
  visitedFaceFlags_.resize(mesh.faces.size());
  visitedFaceIndices_.resize(mesh.faces.size());
  for (unsigned fId = 0; fId < mesh.faces.size(); ++fId) {
    visitedFaceFlags_[fId] = false;
    visitedFaceIndices_[fId] = -1;
  }
  numVisitedFaces_ = 0;

  for (unsigned faceId = 0; faceId < mesh.faces.size(); ++faceId) {
    const TriFace &face = mesh.faces[faceId];
    const Vector3f &v0 = mesh.vertices[face.indices[0]].pos;
    const Vector3f &v1 = mesh.vertices[face.indices[1]].pos;
    const Vector3f &v2 = mesh.vertices[face.indices[2]].pos;

    dis2Leaf_ = maxDepth_;
    lowerCornerOff_X_ = 0;
    lowerCornerOff_Y_ = 0;
    lowerCornerOff_Z_ = 0;
    CreateTree(faceId, v0, v1, v2, root_);
  }
}

void RayMeshIntersection::CreateTree(const unsigned &faceId,
  const Vector3f &v0,
  const Vector3f &v1,
  const Vector3f &v2,
  OctreeNode3D *node) {

  int old_dis2Leaf = dis2Leaf_;
  unsigned old_lowerCornerOff_X = lowerCornerOff_X_;
  unsigned old_lowerCornerOff_Y = lowerCornerOff_Y_;
  unsigned old_lowerCornerOff_Z = lowerCornerOff_Z_;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k  = 0; k < 2; ++k) {
        dis2Leaf_ = old_dis2Leaf - 1;
        unsigned width =  1<<dis2Leaf_;
        lowerCornerOff_X_ = old_lowerCornerOff_X + i*width;
        lowerCornerOff_Y_ = old_lowerCornerOff_Y + j*width;
        lowerCornerOff_Z_ = old_lowerCornerOff_Z + k*width;

        Box3D box;
        box.corners[0][0] = lowerCorner_[0] + lowerCornerOff_X_*leafNodeSize_[0];
        box.corners[0][1] = lowerCorner_[1] + lowerCornerOff_Y_*leafNodeSize_[1];
        box.corners[0][2] = lowerCorner_[2] + lowerCornerOff_Z_*leafNodeSize_[2] ;
        box.corners[1][0] = lowerCorner_[0] + (lowerCornerOff_X_+width)*leafNodeSize_[0];
        box.corners[1][1] = lowerCorner_[1] + (lowerCornerOff_Y_+width)*leafNodeSize_[1];
        box.corners[1][2] = lowerCorner_[2] + (lowerCornerOff_Z_+width)*leafNodeSize_[2];

        if (!TriangleBoxIntersection(v0, v1, v2, box))
          continue;

        unsigned nodeId = 4*i + 2*j + k;

        if (dis2Leaf_ > 0) {
          if (node->m_children[nodeId] == NULL)
            node->m_children[nodeId] = new OctreeNode3D();
          CreateTree(faceId, v0, v1, v2, (OctreeNode3D*)node->m_children[nodeId]);
        } else {
          if (node->m_children[nodeId] == NULL)
            node->m_children[nodeId] = new OctreeLeaf3D();
          

          OctreeLeaf3D *leaf = (OctreeLeaf3D*)(node->m_children[nodeId]);
          leaf->intersectingFaceIndices.push_back(faceId);
        }
      }
    }
  }
}
