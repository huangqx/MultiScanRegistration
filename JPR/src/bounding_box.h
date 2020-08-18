#ifndef bounding_box_h_
#define bounding_box_h_
/**************************************************************************/
// This file implements the data structure for a bounding box of a surface
// represented as a point cloud or a triangular mesh
/**************************************************************************/


#include "linear_algebra.h"

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
  // The point is represented as a float m_vPos[3]
	void Insert_A_Point(float* m_vPos);
  // The point is represented as a 3D vector
	void Insert_A_Point(const Vector3f	&p);
  // The point is represented as a 3D vector
  void Insert_A_Point(const Vector3d& p);
	
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
#endif