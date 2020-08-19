/*
Copyright (c) 2012, Qixing Huang
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Stanford University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef tri_mesh_h_
#define tri_mesh_h_

#include "linear_algebra.h"
#include <vector>
#include <string>
using namespace std;

// Copy the misc.h info here
const int VISIBLE_BIT1= 1<<24;

int Index_(const int &flags);
void SetIndex_(int id, int *flags);
bool Visible_(const int &flags);
void SetVisibleBit_(bool f, int *flags);
Vector3f FaceNormal(const Vector3f &p1, const Vector3f &p2, const Vector3f &p3);

struct BoundingBox {
 public:
	 BoundingBox();
	 ~BoundingBox();

	 /************************************************************************/
	 /* IO Operations
	 */
	 /************************************************************************/
	 void Read(FILE *file_ptr);
	 void Write(FILE *file_ptr);

	 // Update the size of the bounding box 
  // when adding a new bounding box
	 void Insert_A_Box(BoundingBox* newbox);
	
	 // Update the size of the bounding box 
  // when adding a new point
	 void Insert_A_Point(float* m_vPos);

	 void Insert_A_Point(const Vector3f	&p);

  void Insert_A_Point(double x, double y, double z);
	
	 //Insert the feature values of a point
	 void Insert_A_Feature(float* signatures);

	 //Initialize 
	 void Initialize();

	 //The High Corner of this bounding box
	 Vector3d	upper_corner;
	
	 //The Low Corner of this bounding box
	 Vector3d	lower_corner;
	
	 //The size of this bounding box
	 Vector3d	size;

	 //The middle of thie bounding box
	 Vector3d	center_point;
};


struct TriVertex {
 public:
  TriVertex() {
    pos[0] = pos[1] = pos[2] = 0.f;
    nor[0] = nor[1] = nor[2] = 0.f;
    index = 0;
  }
  ~TriVertex() {
  }
  Vector3f pos;
  Vector3f nor;
  int index;
};

struct TriFace {
 public:
  TriFace() {
    indices[0] = indices[1] = indices[2] = 0;
    nor[0] = nor[1] = nor[2] = 0.f;
    flags = 0;
  }
  ~TriFace() {
  }
  void SetVisibleBit(bool f) {
    SetVisibleBit_(f, &flags);
  }
  bool Visible() const {
    return Visible_(flags);
  }
  void SetMeshIndex(int id) {
    SetIndex_(id, &flags);
  }
  int GetMeshIndex() const {
    return Index_(flags);
  }
  int indices[3];
  Vector3f nor;
  int flags;
};

struct TriMesh {
 public:
  TriMesh() {
  }
  ~TriMesh() {
  }
  void Read(FILE *file_ptr);
  void Write(FILE *file_ptr);
  void ComputeBoundingBox();
  void ComputeNormals();
  void Normalize();
  void RemoveUnusedVertices();
  void InitializeFromMatlab(
	  const double *vertexData,
	  const unsigned& numVertices,
	  double *faceData,
	  const unsigned& numFaces);
  vector<TriVertex> vertices;
  vector<TriFace> faces;
  BoundingBox b_box;
};


#endif