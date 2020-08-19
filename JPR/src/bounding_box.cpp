#include "bounding_box.h"
#include "linear_algebra.h"
//#include "stdafx.h"

BoundingBox::BoundingBox() {
	for(int i = 0; i < 3; i++) {
		upper_corner[i] = -1e100;
		lower_corner[i] = 1e100;
    center_point[i] = 0.0;
		size[i] = 0.0;
	}
}

BoundingBox::~BoundingBox() {
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


void BoundingBox::Insert_A_Point(const Vector3d& p) {
  for (int i = 0; i < 3; i++) {
    if (p[i] < lower_corner[i])
      lower_corner[i] = p[i];
    if (p[i] > upper_corner[i])
      upper_corner[i] = p[i];
  }

  center_point = (upper_corner + lower_corner) * 0.5;
  size = upper_corner - lower_corner;
}

void	BoundingBox::Initialize() {
	for(int i = 0; i < 3; i++) {
		upper_corner[i] = -1e100;
		lower_corner[i] = 1e100;
	}
}
