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

#include <algorithm>
using namespace std;

#include "tri_mesh.h"
#include "math.h"

void TriMesh::Normalize() {
  b_box.Initialize();
  for (unsigned fId = 0; fId < faces.size(); ++fId) {
    const TriVertex &v0 = vertices[faces[fId].indices[0]];
    const TriVertex &v1 = vertices[faces[fId].indices[1]];
    const TriVertex &v2 = vertices[faces[fId].indices[2]];
    b_box.Insert_A_Point(v0.pos);
    b_box.Insert_A_Point(v1.pos);
    b_box.Insert_A_Point(v2.pos);
  }
  
  Vector3f center;
  center[0] = static_cast<float> (b_box.center_point[0]);
  center[1] = static_cast<float> (b_box.lower_corner[1]);
  center[2] = static_cast<float> (b_box.center_point[2]);
  float scale = float(b_box.size[1]);

  for (unsigned vId = 0; vId < vertices.size(); ++vId) {
    for (int i = 0; i < 3; ++i) {
      vertices[vId].pos[i] = (vertices[vId].pos[i] - center[i])/scale;
    }
  }
  ComputeBoundingBox();
}

void TriMesh::ComputeNormals() {
  for (unsigned vId = 0; vId < vertices.size(); ++vId) {
    vertices[vId].nor[0] = vertices[vId].nor[1] = vertices[vId].nor[2] = 0.f;
  }

  for (unsigned fId = 0; fId < faces.size(); ++fId) {
    const TriVertex &v0 = vertices[faces[fId].indices[0]];
    const TriVertex &v1 = vertices[faces[fId].indices[1]];
    const TriVertex &v2 = vertices[faces[fId].indices[2]];

    Vector3f nor = FaceNormal(v0.pos, v1.pos, v2.pos);
    faces[fId].nor = nor;
    vertices[faces[fId].indices[0]].nor += nor;
    vertices[faces[fId].indices[1]].nor += nor;
    vertices[faces[fId].indices[2]].nor += nor;
  }

  for (unsigned vId = 0; vId < vertices.size(); ++vId) {
    if (vertices[vId].nor.getSqrNorm() > 1e-10) {
      vertices[vId].nor /= sqrt(vertices[vId].nor*vertices[vId].nor);
    }
  }
}

void TriMesh::ComputeBoundingBox() {
  b_box.Initialize();
  for (unsigned fId = 0; fId < faces.size(); ++fId) {
    const TriVertex &v0 = vertices[faces[fId].indices[0]];
    const TriVertex &v1 = vertices[faces[fId].indices[1]];
    const TriVertex &v2 = vertices[faces[fId].indices[2]];
    b_box.Insert_A_Point(v0.pos);
    b_box.Insert_A_Point(v1.pos);
    b_box.Insert_A_Point(v2.pos);
  }
}

void TriMesh::RemoveUnusedVertices() {
  vector<unsigned> usedFlags;
  usedFlags.resize(vertices.size());
  for (unsigned vId = 0; vId < vertices.size(); ++vId)
    usedFlags[vId] = false;

  for (unsigned fId = 0; fId < faces.size(); ++fId) {
    const TriFace &face = faces[fId];
    for (int i = 0; i < 3; ++i)
      usedFlags[face.indices[i]] = true;
  }

  vector<int> newVIds;
  newVIds.resize(vertices.size());
  for (unsigned vId = 0; vId < vertices.size(); ++vId)
    newVIds[vId] = -1;

  unsigned new_vId = 0;
  for (unsigned vId = 0; vId < vertices.size(); ++vId) {
    if (usedFlags[vId]) {
      vertices[new_vId] = vertices[vId];
      newVIds[vId] = new_vId;
      new_vId++;
    }
  }
  vertices.resize(new_vId);
  for (unsigned fId = 0; fId < faces.size(); ++fId) {
    for (int i = 0; i < 3; ++i) {
      faces[fId].indices[i] = newVIds[faces[fId].indices[i]];
    }
  }
  ComputeBoundingBox();
}

void TriMesh::InitializeFromMatlab(
	const double *vertexData,
	const unsigned& numVertices,
	double *faceData,
	const unsigned& numFaces) {
	vertices.resize(numVertices);
	for (unsigned vId = 0; vId < numVertices; ++vId) {
		for (int i = 0; i < 3; ++i)
			vertices[vId].pos[i] = vertexData[3*vId+i];
	}
	faces.resize(numFaces);
	for (unsigned fId = 0; fId < numFaces; ++fId) {
		for (int i = 0; i < 3; ++i)
			faces[fId].indices[i] = static_cast<int> (faceData[3*fId+i]-0.5);
	}
}

BoundingBox::BoundingBox() {
  for(int i = 0; i < 3; i++) {
    upper_corner[i] = -1e100;
    lower_corner[i] = 1e100;
    size[i] = 0.0;
  }
}

BoundingBox::~BoundingBox() {
}

void BoundingBox::Read(FILE *file_ptr) {
  fread(&lower_corner[0], sizeof(Vector3d), 1, file_ptr);
  fread(&upper_corner[0], sizeof(Vector3d), 1, file_ptr);
  fread(&size[0], sizeof(Vector3d), 1, file_ptr);
  fread(&center_point[0], sizeof(Vector3d), 1, file_ptr);
}

void BoundingBox::Write(FILE *file_ptr) {
  fwrite(&lower_corner[0], sizeof(Vector3d), 1, file_ptr);
  fwrite(&upper_corner[0], sizeof(Vector3d), 1, file_ptr);
  fwrite(&size[0], sizeof(Vector3d), 1, file_ptr);
  fwrite(&center_point[0], sizeof(Vector3d), 1, file_ptr);
}

void BoundingBox::Insert_A_Box(BoundingBox* newbox) {
  for(int i = 0; i < 3; i++) {
    if(newbox->lower_corner[i] < lower_corner[i])
      lower_corner[i] = newbox->lower_corner[i];
    if(newbox->upper_corner[i] > upper_corner[i])
      upper_corner[i] = newbox->upper_corner[i];
  }

  center_point = (upper_corner + lower_corner) * 0.5;
  size = upper_corner - lower_corner;
}

void BoundingBox::Insert_A_Point(float *m_vPos) {
  for(int i = 0; i < 3; i++) {
    if(m_vPos[i] < lower_corner[i])
      lower_corner[i] = m_vPos[i];
    if(m_vPos[i] > upper_corner[i])
      upper_corner[i] = m_vPos[i];
  }

  center_point = (upper_corner + lower_corner) * 0.5;
  size = upper_corner - lower_corner;
}

void BoundingBox::Insert_A_Point(const Vector3f	&p) {
  for(int i = 0; i < 3; i++) {
    if(p[i] < lower_corner[i])
      lower_corner[i] = p[i];
    if(p[i] > upper_corner[i])
      upper_corner[i] = p[i];
  }

  center_point = (upper_corner + lower_corner) * 0.5;
  size = upper_corner - lower_corner;
}

void BoundingBox::Insert_A_Point(double x, double y, double z) {
  if(x < lower_corner[0])
    lower_corner[0] = x;
  if(x > upper_corner[0])
    upper_corner[0] = x;

  if(y < lower_corner[1])
    lower_corner[1] = y;
  if(y > upper_corner[1])
    upper_corner[1] = y;

  if(z < lower_corner[2])
    lower_corner[2] = z;
  if(z > upper_corner[2])
    upper_corner[2] = z;
  center_point = (upper_corner + lower_corner) * 0.5;
  size = upper_corner - lower_corner;
}

void	BoundingBox::Initialize() {
  for(int i = 0; i < 3; i++) {
    upper_corner[i] = -1e100;
    lower_corner[i] = 1e100;
  }
}

int Index_(const int &flags) {
  return flags &((1<<16)-1);
}

void SetIndex_(int id, int *flags) {
  int off = id - Index_(*flags);
  *flags += off;
}

bool Visible_(const int &flags) {
  return (flags & VISIBLE_BIT1) == VISIBLE_BIT1;
}

void SetVisibleBit_(bool f, int *flags) {
  if (f) {
    (*flags) = (*flags) | VISIBLE_BIT1;
  } else {
    (*flags) = (*flags) & (~VISIBLE_BIT1);
  }
}

Vector3f FaceNormal(const Vector3f &p1,
  const Vector3f &p2,
  const Vector3f &p3) {
    Vector3f dev32 = p3 - p2, dev21 = p2 - p1;

    Vector3f normal = dev21.crossProduct(dev32);
    if (normal*normal < 1e-20) {
      normal[0] = normal[1] = normal[2] = 0.f;
      normal[2] = 1.f;
    } else {
      normal /= sqrt(normal*normal);
    }
    return normal;
}