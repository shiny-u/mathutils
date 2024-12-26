//
// Created by frongere on 28/12/17.
//

#ifndef MATHUTILS_FUNCTIONS_H
#define MATHUTILS_FUNCTIONS_H

#include "StdVector.h"
#include "VectorGeneration.h"

namespace mathutils {

  template <class Real>
  std::vector<Real> Sinc(const std::vector<Real>& x) {
      unsigned long n = x.size();
      std::vector<Real> out(n);
      for (unsigned long i=0; i<n; i++) {
          if (x[i] == 0.) {
              out[i] = 1.;
          } else {
              out[i] = sin(x[i]) / x[i];
          }
      }
      return out;
  }

  template <class Real>
  std::vector<Real> Sinc(const Real xmin, const Real xmax, const unsigned int n) {
      auto x = linspace<Real>(xmin, xmax, n);
      return Sinc(x);
  }

  template <class Real>
  std::vector<Real> Heaviside(const std::vector<Real>& x, const Real x_step) {
      unsigned long n = x.size();

      std::vector<Real> out(n);
      Real xval;
      for (unsigned long i=0; i<n; i++) {
          xval = x[i];
          if (xval < x_step) {
              out[i] = 0.;
          } else if (xval > x_step) {
              out[i] = 1.;
          } else {
              out[i] = 0.5;
          }
      }
      return out;
  }

  template <class Real>
  std::vector<Real> Door(const std::vector<Real>& x, const Real x_min, const Real x_max, const Real height=1.) {
      auto out = Heaviside(x, x_min);
      out = Minus(out, Heaviside<Real>(x, x_max));
      Scale<Real>(out, height);
      return out;
  }

  template <class Real>
  std::vector<Real> Dirac(const std::vector<Real>& x, const Real x_impulse, const Real height=1.) {
      unsigned long n = x.size();

      Real dx_2 = 0.5 * (x[1] - x[0]);

      std::vector<Real> out(n);
      for (unsigned long i=0; i<n; i++) {
          if (abs(x[i] - x_impulse) < dx_2) { // FIXME: le critere n'est pas robuste
              out[i] = height;
          } else {
              out[i] = 0.;
          }
      }
      return out;
  }

  template <class Real>
  std::vector<Real> Ramp(const std::vector<Real>& x,
                         const Real x_min, const Real x_max,
                         const Real y_left, const Real y_right) {

      unsigned long n = x.size();

      assert(x_min <= x_max);
      assert(x_min >= x[0] && x_max <= x[n-1]);

      Real a = (y_right - y_left) / (x_max - x_min);
      Real b = y_left - a*x_min;

      std::vector<Real> out(n);
      Real xi;
      for (unsigned long i=0; i<n; i++) {
          xi = x[i];
          if (xi <= x_min) {
              out[i] = y_left;
          } else if (x[i] <= x_max) {
              out[i] = a*xi + b;
          } else {
              out[i] = y_right;
          }
      }

      return out;
  }

}  // end namespace mathutils

#endif //MATHUTILS_FUNCTIONS_H
