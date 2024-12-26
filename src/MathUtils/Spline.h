//
// Created by frongere on 29/12/17.
//

#ifndef MATHUTILS_SPLINE_H
#define MATHUTILS_SPLINE_H

#include "unsupported/Eigen/Splines"

namespace mathutils {

    // TODO: voir pour le passage en 2D... (surfaces)
    template <class Real, int Dim, int Deg>
    class Spline {

    private:
        Real x_min;
        Real x_max;

        // Spline of one-dimensional "points."
        Eigen::Spline<Real, Dim> spline_;


    public:
        Spline(Eigen::Matrix<Real, Eigen::Dynamic, 1> const &x_vec,
               Eigen::Matrix<Real, Eigen::Dynamic, 1> const &y_vec)
                : x_min(x_vec.minCoeff()),
                  x_max(x_vec.maxCoeff()),
                // Spline fitting here. X values are scaled down to [0, 1] for this.
                  spline_(Eigen::SplineFitting<Eigen::Spline<Real, 1>>::Interpolate(
                          y_vec.transpose(),
                          // No more than cubic spline, but accept short vectors.

                          std::min<int>(x_vec.rows() - 1, Deg),
                          scaled_values(x_vec)
                          )
                  )
        { }

        Spline(std::vector<Real>& x_vec, std::vector<Real>& y_vec) {
            // Converting std vectors to eigen vectors
            Eigen::VectorXd xvals = Eigen::Map<Eigen::VectorXd>(x_vec.data(), x_vec.size());
            Eigen::VectorXd yvals = Eigen::Map<Eigen::VectorXd>(y_vec.data(), y_vec.size());

            // Creating a Spline object with standard costructor
            Spline<Real, Dim, Deg> spl(xvals, yvals);

            // Copying Spline attributes into the current instance
            x_min = spl.x_min;
            x_max = spl.x_max;
            spline_ = spl.spline_;
        }

        double operator()(Real x) const {
            // x values need to be scaled down in extraction as well.
            return spline_(scaled_value(x))(0);
        }

    private:
        // Helpers to scale X values down to [0, 1]
        Real scaled_value(Real x) const {
            return (x - x_min) / (x_max - x_min);
        }

        Eigen::RowVectorXd scaled_values(Eigen::VectorXd const &x_vec) const {
            return x_vec.unaryExpr([this](Real x) { return scaled_value(x); }).transpose();
        }


    };


}

#endif //MATHUTILS_SPLINE_H
