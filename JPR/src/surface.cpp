#include "surface.h"

Surface::Surface() {
  flags_	= 0;
  for (int i = 0; i < 3; i++) {
    bounding_box_.center_point[i] = 0.0;
    bounding_box_.size[i] = 0.0;
  }

}

Surface::~Surface() {
}

bool Surface::Render() const {
  return (flags_ & SURFACE_RENDER_BIT) == SURFACE_RENDER_BIT;
}

void Surface::SetRenderBit(bool f) {
  if (f) {
    flags_ = flags_ | SURFACE_RENDER_BIT;
  } else {
    flags_ = flags_ & (~SURFACE_RENDER_BIT);
  }
}

bool Surface::Locked() const {
  return (flags_ & SURFACE_LOCKED_BIT) == SURFACE_LOCKED_BIT;
}

void Surface::SetLockedBit(bool f) {
  if (f) {
    flags_ = flags_ | SURFACE_LOCKED_BIT;
  } else {
    flags_ = flags_ & (~SURFACE_LOCKED_BIT);
  }
}

bool Surface::IsQuad() const {
  return (flags_ & SURFACE_ISQUAD_BIT) == SURFACE_ISQUAD_BIT;
}

void Surface::SetIsQuadBit(bool f) {
  if (f) {
    flags_ = flags_ | SURFACE_ISQUAD_BIT;
  } else {
    flags_ = flags_ & (~SURFACE_ISQUAD_BIT);
  }
}