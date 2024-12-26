//
// Created by pierre-yves on 15/01/20.
//

#ifndef MATHUTILS_INTEGRAND_H
#define MATHUTILS_INTEGRAND_H

namespace mathutils {

  /**
  * Class for handling the function (integrand) to be integrated.
  */
  template<typename T>
  class Integrand {

   public:

    /// This function evalutes the integrand at the point x.
    virtual T Evaluate(const Vector3d<double> &x) const = 0;

  };

}

#endif //MATHUTILS_INTEGRAND_H
