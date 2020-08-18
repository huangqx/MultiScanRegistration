#ifndef affine_transformation_h_
#define affine_transformation_h_

#include "linear_algebra.h"

class Affine2d {   
 private:
	// Translation is the vector
	// followed by a dim x dim matrix
	Vector2d v[3];
 public:
	// ---- constructors
	inline Affine2d();

	inline Vector2d& operator[](const unsigned &index);
	inline const Vector2d& operator[](const unsigned &index) const;

	//  ----  operators +, -, *, /
	inline Affine2d operator+(const Affine2d &op) const;
	inline Affine2d operator-(const Affine2d &op) const;
	inline Affine2d operator-() const;
	inline Affine2d operator*(const double &s) const;

	inline Affine2d operator/(const double &s) const;

	//  ---- operators +=, -=, *=, /=
	inline Affine2d operator+=(const Affine2d &op);
	inline Affine2d operator-=(const Affine2d &op);
	inline Affine2d operator*=(const double &op);
	inline Affine2d operator*=(const Affine2d &op);
	inline Affine2d operator/=(const double &op);
	inline Affine2d operator=(const Affine2d &op);

	inline void Initialize();
	inline void SetZero();
	inline double Det();
	// ---- multiplying with vectors
	inline Vector2d operator*(const Vector2d &v) const;
};

class Affine3d {   
 private:
	// Translation is the vector
	// followed by a dim x dim matrix
	Vector3d v[4];
 public:
	// ---- constructors
	inline Affine3d();
	inline Affine3d(double *velocity);

	inline Vector3d& operator[](const unsigned &index);
	inline const Vector3d& operator[](const unsigned &index) const;

	//  ----  operators +, -, *, /
	inline Affine3d operator+(const Affine3d &op) const;
	inline Affine3d operator-(const Affine3d &op) const;
	inline Affine3d operator-() const;
	inline Affine3d operator*(const double &s) const;
	inline Affine3d operator*(const Affine3d &op) const ;

	inline Affine3d operator/(const double &s) const;

	//  ---- operators +=, -=, *=, /=
	inline Affine3d operator+=(const Affine3d &op);
	inline Affine3d operator-=(const Affine3d &op);
	inline Affine3d operator*=(const double &op);
	inline Affine3d operator*=(const Affine3d &op);
	inline Affine3d operator/=(const double &op);
	inline Affine3d operator=(const Affine3d &op);
	inline Affine3d	Inverse();
	inline Affine3d	Inverse() const;

	inline void Initialize();
	inline void Normalize();
	inline void SetZero();
	inline double Det();
  inline void GetData(Vector12d* vec);
  inline void SetData(const Vector12d& vec);
	// ---- multiplying with vectors
	inline Vector3d operator*(const Vector3d &v) const;
};

#endif
