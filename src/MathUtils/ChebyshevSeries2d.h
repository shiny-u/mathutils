//
// Created by pierre-yves on 30/11/2020.
//

#ifndef MATHUTILS_CHEBYSHEVSERIES2D_H
#define MATHUTILS_CHEBYSHEVSERIES2D_H

#include <memory>
#include "Matrix.h"
#include "BoostFunctions.h"

namespace mathutils {

/**
* Class for handling 2d functions for the double Chebyshev series approximation.
*/
  template<typename T>
  class Function2d {

   public:

    /// This function evaluates the function at the point (x, y).
    virtual T Evaluate(const double &x, const double &y) const = 0;

  };

  /**
 * Class for handling and computing the double Chebyshev series approximation of a function.
 */
  template<typename T>
  class ChebyshevSeries2dBase {

   public:

    /// Contructor of the class.
    ChebyshevSeries2dBase(Function2d<T> *F, const double &xmin, const double &ymin, const int &order_x,
                                    const int &order_y) : m_function(F), m_x_min(xmin), m_y_min(ymin),
                                    m_order_x(order_x), m_order_y(order_y) {
      m_aij = MatrixMN<T>::Zero(order_x + 1, order_y + 1);
    }

   public:

    /// This method computes the coefficients aij.
    void Compute_aij() {

      // x abscissa.
      std::vector<double> x_tilde; // In [-1,1].
      std::vector<double> x; // In [xmin,xmax].

      // x.
      for (int r = 0; r <= m_order_x; ++r) {
        double tmp = cos(MU_PI_2 * (2. * r + 1.) / (m_order_x + 1.));
        x_tilde.push_back(tmp);
        x.push_back(AffineTransformationUnitToSegment_x(tmp));
      }

      // y abscissa.
      std::vector<double> y_tilde; // In [-1,1].
      std::vector<double> y; // In [ymin,ymax].

      // y.
      for (int s = 0; s <= m_order_y; ++s) {
        double tmp = cos(MU_PI_2 * (2. * s + 1.) / (m_order_y + 1.));
        y_tilde.push_back(tmp);
        y.push_back(AffineTransformationUnitToSegment_y(tmp));
      }

      // aij.
      // The Clenshaw algorithm cannot be used as x_tilde and y_tilde depend on r and s.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          for (int r = 0; r <= m_order_x; ++r) {
            double Ti = Chebyshev_polynomial(i, x_tilde.at(r));
            for (int s = 0; s <= m_order_y; ++s) {
              m_aij(i, j) += this->m_function->Evaluate(x.at(r), y.at(s)) * Ti * Chebyshev_polynomial(j, y_tilde.at(s));
            }
          }
          m_aij(i, j) *= 4. / ((m_order_x + 1.) * (m_order_y + 1.));
          if (i == 0 and j == 0) {
            m_aij(i, j) /= 4.;
          } else if ((i == 0 and j != 0) or (i != 0 and j == 0)) {
            m_aij(i, j) /= 2.;
          }
        }
      }
    }

    /// This method computes the double Chebyshev series approximation.
    T Evaluate(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      // Partial sum using Clenshaw algorithm.
      std::vector<double> ai;
      ai.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_aij.row(i);
        std::vector<double> bj(tmp.data(), tmp.data() + tmp.size());
        bj.at(0) *= 2.;
        ai.push_back(boost::math::chebyshev_clenshaw_recurrence(bj.data(), bj.size(), yunit));
      }
      ai.at(0) *= 2.;
      T result = boost::math::chebyshev_clenshaw_recurrence(ai.data(), ai.size(), xunit);

      return result;

    }

    /// This method computes the x-derivative double Chebyshev series approximation.
    T Evaluate_derivative_x(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      // TODO: Implementation of the Clenshaw algorithm for derivative to study.

      // Partial sum using Clenshaw algorithm except for the derivative.
      T result = 0.;
      for (int i = 0; i <= m_order_x; ++i) {
        double Ti = Chebyshev_polynomial_derivative(i, xunit);
        Eigen::VectorXd tmp = m_aij.row(i);
        std::vector<double> bj(tmp.data(), tmp.data() + tmp.size());
        bj.at(0) *= 2.;
        result += Ti * boost::math::chebyshev_clenshaw_recurrence(bj.data(), bj.size(), yunit);
      }
      result *= CoefficientDerivative_x(x); // With a closed segment, this coefficient is independent of normal_x.

      return result;

    }

    /// This method computes the y-derivative double Chebyshev series approximation.
    T Evaluate_derivative_y(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      // TODO: Implementation of the Clenshaw algorithm for derivative to study.

      // Partial sum using Clenshaw algorithm except for the derivative.
      T result = 0.;
      for (int j = 0; j <= m_order_y; ++j) {
        double Tj = Chebyshev_polynomial_derivative(j, yunit);
        Eigen::VectorXd tmp = m_aij.col(j);
        std::vector<double> ai(tmp.data(), tmp.data() + tmp.size());
        ai.at(0) *= 2.;
        result += Tj * boost::math::chebyshev_clenshaw_recurrence(ai.data(), ai.size(), xunit);
      }
      result *= CoefficientDerivative_y(y); // With a closed segment, this coefficient is independent of normal_y.

      return result;

    }

    /// This method computes the bij of the power series from the aij of the Chebyshev series.
    MatrixMN<T> Compute_bij() {

      if(m_order_x > 18 or m_order_y > 18) {
        std::cout << "The conversion Chebyshev to power series have large numerical inaccuracies for important order." << std::endl;
        std::cout << "The maximum order in x and y is 18." << std::endl;
        std::cout << "Order in x: " << m_order_x << std::endl;
        std::cout << "Order in y: " << m_order_y << std::endl;
        exit(0);
      }

      // Initialization.
      MatrixMN<T> bij = MatrixMN<T>::Zero(m_order_x + 1, m_order_y + 1);

      // Computation.
      for (int i = 0; i <= m_order_x; ++i) {
        for (int j = 0; j <= m_order_y; ++j) {
          for (int r = i; r <= m_order_x; ++r) {
            double lambda_ri = TransformationCoefficient(r, i);
            double tmp = 0.;
            for (int s = j; s <= m_order_y; ++s) {
              double lambda_sj = TransformationCoefficient(s, j);
              tmp += m_aij(r, s) * lambda_sj;
            }
            bij(i, j) += tmp * lambda_ri;
          }
        }
      }

      return bij;

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

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    virtual double AffineTransformationUnitToSegment_y(const double& yunit) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double& xdomain) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    virtual double AffineTransformationSegmentToUnit_y(const double& ydomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double& xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    virtual double CoefficientDerivative_y(const double& ydomain) const = 0;

   protected:

    /// 2d function to be approximated.
    Function2d<T> *m_function;

    /// Order of the series approximation for x.
    int m_order_x;

    /// Order of the series approximation for y.
    int m_order_y;

    /// x minimum.
    double m_x_min;

    /// y minimum.
    double m_y_min;

    /// aij coefficients.
    MatrixMN<T> m_aij;

  };

  /**
 * Class for computing the double Chebyshev series approximation of a function over [xmin, xmax] x [ymin, ymax].
 */
  template<typename T>
  class ChebyshevSeries2dClosed : public ChebyshevSeries2dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevSeries2dClosed(Function2d<T> *F, const double &xmin, const double &xmax, const double &ymin,
                            const double &ymax, const int &order_x, const int &order_y) :
        ChebyshevSeries2dBase<T>(F, xmin, ymin, order_x, order_y), m_x_max(xmax), m_y_max(ymax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double &yunit) const override {
      return 0.5 * (m_y_max - this->m_y_min) * yunit + 0.5 * (m_y_max + this->m_y_min);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return (2. / (m_y_max - this->m_y_min)) * (ydomain - 0.5 * (m_y_max + this->m_y_min));
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const {
      return 2. / (m_y_max - this->m_y_min);
    }

   private:

    /// x maximum.
    double m_x_max;

    /// y maximum.
    double m_y_max;

  };

  /**
 * Class for computing the double Chebyshev series approximation of a function over [xmin, +infinity] x [ymin, +infinity].
 */
  template<typename T>
  class ChebyshevSeries2dOpened : public ChebyshevSeries2dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevSeries2dOpened(Function2d<T> *F, const double &xmin, const double &ymin, const int &order_x,
                            const int &order_y) :
        ChebyshevSeries2dBase<T>(F, xmin, ymin, order_x, order_y) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 2. * (this->m_x_min / (1 - xunit));
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double &yunit) const override {
      return 2. * (this->m_y_min / (1 - yunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return 1. - 2. * (this->m_x_min / xdomain);
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. * this->m_x_min / (xdomain * xdomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

  };

  /**
 * Class for computing the double Chebyshev series approximation of a function over [xmin, xmax] x [ymin, +infinity].
 */
  template<typename T>
  class ChebyshevSeries2dMixed : public ChebyshevSeries2dBase<T> {

   public:

    /// Contructor of the class.
    ChebyshevSeries2dMixed(Function2d<T> *F, const double &xmin, const double &xmax, const double &ymin,
                           const int &order_x, const int &order_y) :
        ChebyshevSeries2dBase<T>(F, xmin, ymin, order_x, order_y), m_x_max(xmax) {
    }

   private:

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for x.
    double AffineTransformationUnitToSegment_x(const double &xunit) const override {
      return 0.5 * (m_x_max - this->m_x_min) * xunit + 0.5 * (m_x_max + this->m_x_min);
    }

    /// This method applied an affine transformation from [-1, 1] to the domain of the approximation for y.
    double AffineTransformationUnitToSegment_y(const double &yunit) const override {
      return 2. * (this->m_y_min / (1 - yunit));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    double AffineTransformationSegmentToUnit_x(const double &xdomain) const override {
      return (2. / (m_x_max - this->m_x_min)) * (xdomain - 0.5 * (m_x_max + this->m_x_min));
    }

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    double AffineTransformationSegmentToUnit_y(const double &ydomain) const override {
      return 1. - 2. * (this->m_y_min / ydomain);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    double CoefficientDerivative_x(const double &xdomain) const {
      return 2. / (m_x_max - this->m_x_min);
    }

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    double CoefficientDerivative_y(const double &ydomain) const {
      return 2. * this->m_y_min / (ydomain * ydomain);
    }

    /// x maximum.
    double m_x_max;

  };

}

#endif //MATHUTILS_CHEBYSHEVSERIES2D_H
