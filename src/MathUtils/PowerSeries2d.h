//
// Created by pierre-yves on 03/03/2021.
//

#ifndef MATHUTILS_POWERSERIES2D_H
#define MATHUTILS_POWERSERIES2D_H

#include "Matrix.h"
#include "Polynomial.h"
#include "Functions.h"

namespace mathutils {

/**
* Class for handling and computing the double power series approximation.
*/
  template<typename T>
  class PowerSeries2dBase {

   public:

    /// Contructor of the class.
    PowerSeries2dBase(const MatrixMN <T> &bij, const double &xmin, const double &ymin) : m_bij(bij),
                                 m_x_min(xmin), m_y_min(ymin) {
      m_order_x = bij.rows() - 1;
      m_order_y = bij.cols() - 1;
    }

    /// This method computes the double power series approximation for unit coordinates.
    T Evaluate_unit(const double &xunit, const double &yunit) const {

      // Partial sums.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_bij.row(i);
        std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
        qi.push_back(Horner<double>(pj, yunit));
      }
      T result = Horner<double>(qi, xunit);

      return result;

    }

    /// This method computes the double power series approximation.
    T Evaluate(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      return Evaluate_unit(xunit, yunit);

    }

    /// This method computes the x-derivative double power series approximation for unit coordinates.
    T Evaluate_derivative_x_unit(const double &x, const double &xunit, const double &yunit) const {

      // Partial sums.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_bij.row(i);
        std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
        qi.push_back(Horner<double>(pj, yunit));
      }
      T result = CoefficientDerivative_x(x) * Horner_derivative<double>(qi, xunit);

      return result;

    }

    /// This method computes the x-derivative double power series approximation.
    T Evaluate_derivative_x(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      return Evaluate_derivative_x_unit(x, xunit, yunit);

    }

    /// This method computes the y-derivative double power series approximation for unit coordinates.
    T Evaluate_derivative_y_unit(const double &y, const double &xunit, const double &yunit) const {

      // Partial sums.
      std::vector<double> qi;
      qi.reserve(m_order_x + 1);
      for (int i = 0; i <= m_order_x; ++i) {
        Eigen::VectorXd tmp = m_bij.row(i);
        std::vector<double> pj(tmp.data(), tmp.data() + tmp.size());
        qi.push_back(Horner_derivative<double>(pj, yunit));
      }
      T result = CoefficientDerivative_y(y) * Horner<double>(qi, xunit);

      return result;

    }

    /// This method computes the y-derivative double power series approximation.
    T Evaluate_derivative_y(const double &x, const double &y) const {

      // Parameters.
      double xunit = AffineTransformationSegmentToUnit_x(x);
      double yunit = AffineTransformationSegmentToUnit_y(y);

      return Evaluate_derivative_y_unit(y, xunit, yunit);

    }

   protected:

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for x.
    virtual double AffineTransformationSegmentToUnit_x(const double &xdomain) const = 0;

    /// This method applied an affine transformation from the domain of the approximation to [-1, 1] for y.
    virtual double AffineTransformationSegmentToUnit_y(const double &ydomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt x.
    virtual double CoefficientDerivative_x(const double &xdomain) const = 0;

    /// This method gives the coefficient in front of the derivate of the  double Chebychev series wrt y.
    virtual double CoefficientDerivative_y(const double &ydomain) const = 0;

   protected:

    /// Order of the series approximation for x.
    int m_order_x;

    /// Order of the series approximation for y.
    int m_order_y;

    /// x minimum.
    double m_x_min;

    /// y minimum.
    double m_y_min;

    /// Pij coefficients.
    MatrixMN<T> m_bij;

  };

  /**
  * Class for computing the double power series approximation over [xmin, xmax] x [ymin, ymax].
  */
  template<typename T>
  class PowerSeries2dClosed : public PowerSeries2dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries2dClosed(const MatrixMN <T> &bij, const double &xmin, const double &xmax, const double &ymin,
                        const double &ymax) : PowerSeries2dBase<T>(bij, xmin, ymin), m_x_max(xmax), m_y_max(ymax) {}

   private:

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
  * Class for computing the double power series approximation of a function over [xmin, +infinity] x [ymin, +infinity].
  */
  template<typename T>
  class PowerSeries2dOpened : public PowerSeries2dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries2dOpened(const MatrixMN <T> &bij, const double &xmin, const double &ymin) :
    PowerSeries2dBase <T>(bij, xmin, ymin) {}

   private:

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
  * Class for computing the double power series approximation of a function over [xmin, xmax] x [ymin, +infinity].
  */
  template<typename T>
  class PowerSeries2dMixed : public PowerSeries2dBase<T> {

   public:

    /// Contructor of the class.
    PowerSeries2dMixed(const MatrixMN <T> &bij, const double &xmin, const double &xmax, const double &ymin) :
    PowerSeries2dBase<T>(bij, xmin, ymin), m_x_max(xmax) {}

   private:

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

#endif //MATHUTILS_POWERSERIES2D_H