#ifndef latent_surf_h_
#define latent_surf_h_

#include "linear_algebra.h"
#include "affine_transformation.h"
#include "bounding_box.h"
#include <vector>
using namespace std;


class LatentSurf {
 public:
	LatentSurf() {
  }
	~LatentSurf() {
  }

	/************************************************************************/
	/* IO Operations
	*/
	/************************************************************************/
	void Read(FILE *file_ptr);
	void Write(FILE *file_ptr);

  BoundingBox b_box;
};
#endif