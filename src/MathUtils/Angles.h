//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_ANGLES_H
#define MATHUTILS_ANGLES_H

#include <string>
#include <iostream>

#include "Unit.h"
#include "Constants.h"

namespace mathutils {


  /// Getting the trigonometric quantities (cosine and sine) from an angle in one pass with unit management
  template<class Scalar=double>
  inline void _GetTrigo(const Scalar pangle, Scalar &cosAngle, Scalar &sinAngle, ANGLE_UNIT unit = RAD) {
    Scalar angle = pangle;
    if (unit == DEG) {
      angle *= MU_PI_180;
    }
    cosAngle = cos(angle);
    sinAngle = sin(angle);
  }

  /// Reminder of mod(2*pi) to put back an angle expressed in radian into [0, 2pi[
  template<class Real=double>
  inline Real modulo2pi(const Real a) {
    return fmod(a, (Real) MU_2PI);
  }

  /// Reminder of mod(360) to put back an angle expressed in degrees into [0, 360[
  template<class Real=double>
  inline Real modulo360(const Real a) {
    return fmod(a, (Real) 360.);
  }

  /// Normalizes angles to [0, 360[ range
  template<class Real>
  inline Real Normalize_0_360(const Real a) {
    Real angle = modulo360(a);
    if (angle < 0) {
      angle += 360;
    }
    return angle;
  }

  /// Normalizes angles to ]-180, 180]
  template<class Real>
  inline Real Normalize__180_180(const Real a) {
    Real angle = modulo360(a + 180);
    if (angle <= 0) {
      angle += 360;
    }
    return angle - 180;
  }

  /// Normalizes angles to [0, 2*PI[ range
  template<class Real>
  inline Real Normalize_0_2PI(const Real a) {
    Real angle = modulo2pi(a);
    if (angle < 0) {
      angle += MU_2PI;
    }
    return angle;
  }

  /// Normalizes angles to ]-PI, PI]
  template<class Real>
  inline Real Normalize__PI_PI(const Real a) {
    Real angle = modulo2pi(a + M_PI);
    if (angle <= 0) {
      angle += MU_2PI;
    }
    return angle - M_PI;
  }

  /// Convert string to angular unit
  inline ANGLE_UNIT STRING2ANGLE(const std::string unit) {
    ANGLE_UNIT angle_unit;
    if (unit == "DEG") {
      angle_unit = DEG;
    } else if (unit == "RAD") {
      angle_unit = RAD;
    } else {
      std::cerr << "unknown value for the unit angle" << std::endl;
      exit(1);
    }
    return angle_unit;
  }


}

#endif //MATHUTILS_UTILS_H
