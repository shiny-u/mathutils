//
// Created by pierre-yves on 30/11/2020.
//

#ifndef MATHUTILS_CHEBYSHEVSERIES1D_H
#define MATHUTILS_CHEBYSHEVSERIES1D_H

#include <memory>
#include "Matrix.h"
#include "BoostFunctions.h"

namespace mathutils {

  /**
  * Class for handling 1d functions for the single Chebyshev series approximation.
  */
  template<typename T>
  class Function1d {

   public:

    /// This function evaluates the function at the point x.
    virtual T Evaluate(const double &x) const = 0;

  };

  /**
 * Class for handling and computing the single Chebyshev series approximation of a function.
 */
  template<typename T>
  class ChebyshevSeries1dBase {

   public:

    /// Contructor of the class.
    ChebyshevSeries1dBase(Function1d<T> *F, const double &xmin, const int &order_x) : m_function(F), m_x_min(xmin),
                          m_order_x(order_x) {
      m_ai = VectorN<T>::Zero(order_x + 1);
    }

   public:

    /// This method computes the coefficients ai.
    void Compute_ai() {

      // x abscissa.
      std::vector<double> x_tilde; // In [-1,1].
      std::vector<double> x; // In [xmin,xmax].

      // x.
      for (int r = 0; r <= m_order_x; ++r) {
        double tmp = cos(MU_PI_2 * (2. * r + 1.) / (m_order_x + 1.));
        x_tilde.push_back(tmp);
        x.push_back(AffineTransformationUnitToSegment_x(tmp));
      }

      // aij.
      // The Clenshaw algorithm cannot be used as x_tilde and y_tilde depend on r and s.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int r = 0; r <= m_order_x; ++r) {
          double Ti = Chebyshev_polynomial(i, x_tilde.at(r));
          m_ai(i) += this->m_function->Evaluate(x.at(r)) * Ti;
        }
        m_ai(i) *= 2. / (m_order_x + 1.);
        if (i == 0) {
          m_ai(i) /= 2.;
        }
      }
    }

    /// This method computes the single Chebyshev series approximation.
    T Evaluate(const double &x) const {

      // Parameter.
      double xunit = AffineTransformationSegmentToUnit_x(x);

      // Partial sum using Clenshaw algorithm.
      std::vector<double> ai;
      ai.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        ai.push_back(m_ai(i));
      }
      ai.at(0) *= 2.;
      T result = boost::math::chebyshev_clenshaw_recurrence(ai.data(), ai.size(), xunit);

      return result;

    }

    /// This method computes the x-derivative single Chebyshev series approximation.
    T Evaluate_derivative_x(const double &x) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);

      // TODO: Implementation of the Clenshaw algorithm for derivative to study.

      // Partial sum using Clenshaw algorithm except for the derivative.
      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial_derivative(i, xunit);
        result += Ti * m_ai(i);
      }
      result *= CoefficientDerivative_x(x); // With a closed segment, this coefficient is independent of normal_x.

      return result;

    }

    /// This method computes the bi of the power series from the ai of the Chebyshev series.
    VectorN<T> Compute_bi() {

      if(m_order_x > 18) {
        std::cout << "The conversion Chebyshev to power series have large numerical inaccuracies for important order." << std::endl;
        std::cout << "The maximum order in xis 18." << std::endl;
        std::cout << "Order in x: " << m_order_x << std::endl;
        exit(0);
      }

      // Initialization.
      VectorN<T> bi = VectorN<T>::Zero(m_order_x + 1);

      // Computation.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int r = i; r <= m_order_x; ++r) {
          double lambda_ri = TransformationCoefficient(r, i);
          bi(i) += m_ai(r) * lambda_ri;
        }
      }

      return bi;

    }

   private:

    /// This method computes the transformation coefficient.
    double TransformationCoefficient(const int& lj, const int& kj) {
      double lambda;
      int ljkj_plus = lj + kj;
      if (lj == 0 and kj == 0) {
        lambda = 1.;
      } else if (ljkj_plus % 2 == 1) { // (lj + kj) impair.
        lambda = 0.;
      } else { // (lj + kj) pair.
        int floor_plus = floor(0.5 * ljkj_plus);
        int floor_minus = floor(0.5 * (lj - kj));
        lambda = pow(-1, floor_minus) * pow(2., kj - 1.) * lj * (Factorial<int, double>(floor_plus - 1) /
                 (Factorial<int, double>(floor_minus) * Factorial<int, double>(kj)));
      }

      return lambda;

    }

   protected:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    virtual double AffineTransformationUnitToSegment_x(const double& xunit) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double& xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the single Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double& xdomain) const = 0;

   protected:

    /// 2d function to be approximated.
    Function1d<T> *m_function;

    /// Order of the series approximation for x.
    int m_order_x;

    /// x minimum.
    double m_x_min;

    /// aij coefficients.
    VectorN<T> m_ai;

  };

  /**
 * Class for computing the single Chebyshev series approximation of a function over [xmin, xmax].
 */
  template<typename T>
  class ChebyshevSeries1dClosed : public ChebyshevSeries1dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevSeries1dClosed(Function1d<T> *F, const double &xmin, const double &xmax, const int &order_x) :
                            ChebyshevSeries1dBase<T>(F, xmin, order_x), m_x_max(xmax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

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
 * Class for computing the single Chebyshev series approximation of a function over [xmin, +infinity].
 */
  template<typename T>
  class ChebyshevSeries1dOpened : public ChebyshevSeries1dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevSeries1dOpened(Function1d<T> *F, const double &xmin, const int &order_x) :
                            ChebyshevSeries1dBase<T>(F, xmin, order_x) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 2. * (this->m_x_min / (1 - xunit));
    }

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

#endif //MATHUTILS_CHEBYSHEVSERIES1D_H
