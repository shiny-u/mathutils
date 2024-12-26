//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2D_H
#define MATHUTILS_INTEGRATION2D_H

#include "QuadratureTables.h"
#include <memory>
#include "Vector3d.h"
#include "Integrand.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations based on a Gauss quadrature.
  */
  template<typename T>
  class Integration2d {

   public:

    /// Contructor of the class.
    Integration2d(Integrand<T>* F, const int& order) : m_order(order), m_integrand(F) {
      m_tables = std::make_unique<QuadratureTables>();
    }

    /// This function computes the 2d integration.
    virtual T Compute(const Vector3d<double> &vertex_1, const Vector3d<double> &vertex_2,
                         const Vector3d<double> &vertex_3) const = 0;

   protected:

    /// Function to be integrated
    Integrand<T>* m_integrand;

    /// Order of the quadrature.
    int m_order;

    /// Table of quadrature coefficients.
    std::unique_ptr<QuadratureTables> m_tables;

  };

}

#endif //MATHUTILS_INTEGRATION2D_H
