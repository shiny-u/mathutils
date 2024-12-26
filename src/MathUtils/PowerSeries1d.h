//
// Created by pierre-yves on 03/03/2021.
//

#ifndef MATHUTILS_POWERSERIES1D_H
#define MATHUTILS_POWERSERIES1D_H

#include "Matrix.h"
#include "Polynomial.h"
#include "Functions.h"

namespace mathutils {

/**
* Class for handling and computing the single power series approximation.
*/
  template<typename T>
  class PowerSeries1dBase {

   public:

    /// Contructor of the class.
    PowerSeries1dBase(const VectorN<T> &bi, const double &xmin) : m_bi(bi), m_x_min(xmin) {
      m_order_x = bi.rows() - 1;
    }

    /// This method computes the single power series approximation for unit coordinates.
    T Evaluate_unit(const double &xunit) const {

      // Partial sum.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        qi.push_back(m_bi(i));
      }
      T result = Horner<double>(qi, xunit);

      return result;

    }

    /// This method computes the single power series approximation.
    T Evaluate(const double &x) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);

      return Evaluate_unit(xunit);

    }

    /// This method computes the x-derivative single power series approximation for unit coordinates.
    T Evaluate_derivative_x_unit(const double &x, const double &xunit) const {

      // Partial sum.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        qi.push_back(m_bi(i));
      }
      T result = CoefficientDerivative_x(x) * Horner_derivative<double>(qi, xunit);

      return result;

    }

    /// This method computes the x-derivative single power series approximation.
    T Evaluate_derivative_x(const double &x) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);

      return Evaluate_derivative_x_unit(x, xunit);

    }

   protected:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double &xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double &xdomain) const = 0;

   protected:

    /// Order of the series approximation for x.
    int m_order_x;

    /// x minimum.
    double m_x_min;

    /// Pij coefficients.
    VectorN<T> m_bi;

  };

  /**
  * Class for computing the single power series approximation over [xmin, xmax].
  */
  template<typename T>
  class PowerSeries1dClosed : public PowerSeries1dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries1dClosed(const VectorN<T> &bi, const double &xmin, const double &xmax) :
                        PowerSeries1dBase<T>(bi, xmin), m_x_max(xmax) {}

   private:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method gives the coefficient in front of the derivate of the single Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. / (m_x_max - this->m_x_min);
    }

   private:

    /// x maximum.
    double m_x_max;

  };

  /**
  * Class for computing the single power series approximation of a function over [xmin, +infinity].
  */
  template<typename T>
  class PowerSeries1dOpened : public PowerSeries1dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries1dOpened(const VectorN<T> &bi, const double &xmin) : PowerSeries1dBase<T>(bi, xmin) {}

   private:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return 1. - 2. * (this->m_x_min / xdomain);
    }

    /// This method gives the coefficient in front of the derivate of the single Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. * this->m_x_min / (xdomain * xdomain);
    }

  };

}

#endif //MATHUTILS_POWERSERIES1D_H