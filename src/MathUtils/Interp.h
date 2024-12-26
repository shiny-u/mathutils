//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_INTERP_H
#define MATHUTILS_INTERP_H

#include <vector>
#include <memory>
#include <cassert>
#include <algorithm>

namespace mathutils {

  enum INTERP_METHOD {
    LINEAR,
    LINEAR_EXTRAPOLATE,
    LINEAR_SATURATE,
    BSPLINE
  };

}  // end namespace mathutils


#endif //MATHUTILS_INTERP_H
