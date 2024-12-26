//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_UNIT_H
#define MATHUTILS_UNIT_H

// TODO: voir à intégrer le framework de gestion des unités :   github.com/nholthaus/units

#include <iostream>  // TODO: retirer et mettre en place un systeme d'erreur partout ou c'est utilise...
#include <vector>
#include "Constants.h"


/*
 * Here we define unit related functions such as conversions between different units
 */

namespace mathutils {


  /* ===========================================================================================
   * ANGLES Conversions functions
   */

  enum ANGLE_UNIT {
    DEG,
    RAD
  };

  /// Conversion DEG->RAD inplace
  template<class T>
  inline void deg2rad(std::vector<T> &vector) {
    for (auto &v : vector) v *= DEG2RAD;
  }

  /// CONVERSION DEG->RAD by copy
  template<class T>
  inline void deg2rad(const std::vector<T> &vector, std::vector<T> &out) {
    out = vector;
    return deg2rad(out);
  }

  /// Conversion RAD->DEG inplace
  template<class T>
  inline void rad2deg(std::vector<T> &vector) {
    for (auto &v : vector) v *= RAD2DEG;
  }

  /// CONVERSION RAD->DEG by copy
  template<class T>
  // TODO : faire la meme inplace...
  inline std::vector<T> rad2deg(const std::vector<T> &vector, std::vector<T> &out) {
    out = vector;
    rad2deg(vector);
  }



//    /// CONVERSION RAD->DEG
//    template <class T>
//    inline T degrees(const T a) {
//        return a * MU_180_PI;
//    }
//
//    // TODO: both two following functions should be realized as MACROS for performance
//    /// CONVERSION DEG->RAD
//    template <class Real=double>
//    inline Real deg2rad(const Real a) {
//        return a * MU_PI_180;
//    }
//
//    /// CONVERSION RAD->DEG
//    template <class Real=double>
//    inline Real rad2deg(const Real a) {
//        return a * MU_180_PI;
//    }

  /* ===========================================================================================
   * Velocities conversions functions
   */

  enum SPEED_UNIT {
    MS,     ///> M/S
    KNOT,   ///> NAUTICAL KNOTS
    KMH     ///> KM/H
  };


  /// Convert nautical knots into m/s
  template<class Real=double>
  inline Real KNOT2MS(Real velocity) {
    return velocity * MU_KNOT;
  }

  /// Convert m/s into nautical knots
  template<class Real=double>
  inline Real MS2KNOT(Real velocity) {
    return velocity / MU_KNOT;
  }

  /// Convert km/h into m/s
  template<class Real=double>
  inline Real KMH2MS(Real velocity) {
    return velocity * MU_KMH;
  }

  /// Convert m/s into km/h
  template<class Real=double>
  inline Real MS2KMH(Real velocity) {
    return velocity / MU_KMH;
  }

  /// Convert km/h into nautical knots
  template<class Real=double>
  inline Real KMH2KNOT(Real velocity) {
    return MS2KNOT(KMH2MS(velocity));
  }

  /// Convert nautical knots into km/h
  template<class Real=double>
  inline Real KNOT2KMH(Real velocity) {
    return MS2KMH(KNOT2MS(velocity));
  }

  /// Conversion kt2ms inplace
  template<class T>
  inline void kt2ms(std::vector<T> &vector) {
    for (auto &v : vector) v *= MU_KNOT;
  }

  /// CONVERSION kt2ms by copy
  template<class T>
  inline void kt2ms(const std::vector<T> &vector, std::vector<T> &out) {
    out = vector;
    return kt2ms(out);
  }

  /// Converts a velocity from a unit to another
  template<class Real=double>
  inline Real convert_velocity_unit(const Real velocity, SPEED_UNIT current_unit, SPEED_UNIT new_unit = MS) {

    Real new_vel;
    // Expressing in M/S
    switch (current_unit) {
      case MS:
        new_vel = velocity;
        break;
      case KMH:
        new_vel = KMH2MS(velocity);
        break;
      case KNOT:
        new_vel = KNOT2MS(velocity);
        break;
    }

    // EXPRESSING IN NEW UNIT
    Real out;
    switch (new_unit) {
      case MS:
        out = new_vel;
        break;
      case KMH:
        out = MS2KMH(new_vel);
        break;
      case KNOT:
        out = MS2KNOT(new_vel);
        break;
    }
    return out;

  }


  /* ===========================================================================================
   * Frequency conversion functions
   */

  enum FREQUENCY_UNIT {
    HERTZ,  ///< Hertz (s**-1)
    RADS,   ///< rad/s
    DEGS,   ///< deg/s
    DEGM,   ///< deg/min
    S,      ///< seconds (period)
    RPM     ///< Round Per Minute
  };  // TODO: ajouter la conversion avec RPM ...

  template<class Real=double>
  inline Real DEGS2RADS(const Real degs) {
    return degs * MU_PI_180;
  }

  template<class Real=double>
  inline Real RADS2DEGS(const Real rads) {
    return rads * MU_180_PI;
  }


  template<class Real=double>
  inline Real HZ2RADS(const Real hz) {
    return MU_2PI * hz;
  }

  template<class Real=double>
  inline Real RADS2HZ(const Real rads) {
    return rads / MU_2PI;
  }

  template<class Real=double>
  inline Real HZ2S(const Real hz) {
    return 1. / hz;
  }

  template<class Real=double>
  inline Real S2HZ(const Real s) {
    return 1. / s;
  }

  template<class Real=double>
  inline Real RADS2S(const Real rads) {
    return MU_2PI / rads;
  }

  template<class Real=double>
  inline Real S2RADS(const Real s) {
    return MU_2PI / s;
  }

  template<class Real=double>
  inline Real RPM2RADS(const Real rpm) {
    return rpm * MU_RPM2RADS;
  }

  template<class Real=double>
  inline Real RADS2RPM(const Real rads) {
    return rads / MU_RPM2RADS;
  }

  template<class Real=double>
  inline Real DEGM2RADS(const Real degm) {
    return degm * MU_PI_180 / 60;
  }

  template<class Real=double>
  inline Real RADS2DEGM(const Real rads) {
    return rads * MU_180_PI * 60;
  }


  template<class Real=double>
  inline Real convert_frequency(const Real in, FREQUENCY_UNIT src_unit, FREQUENCY_UNIT target_unit) {

    if (src_unit == target_unit) return in;

    Real piv;  // Must be in Hz

    switch (src_unit) {
      case HERTZ:
        piv = in;
        break;
      case RADS:
        piv = RADS2HZ(in);
        break;
      case S:
        piv = S2HZ(in);
        break;
      case DEGS:
        piv = RADS2HZ(DEGS2RADS(in));
        break;
      case DEGM:
        piv = RADS2HZ(DEGM2RADS(in));
        break;
      case RPM:
        piv = RADS2HZ(RPM2RADS(in));
        break;
    }

    Real out;
    switch (target_unit) {
      case RADS:
        out = HZ2RADS(piv);
        break;
      case S:
        out = HZ2S(piv);
        break;
      case HERTZ:
        out = piv;
        break;
      case DEGS:
        out = RADS2DEGS(HZ2RADS(piv));
        break;
      case DEGM:
        out = RADS2DEGM(HZ2RADS(piv));
        break;
      case RPM:
        out = RADS2RPM(HZ2RADS(piv));
        break;
    }
    return out;
  }


}

#endif //MATHUTILS_UNIT_H
