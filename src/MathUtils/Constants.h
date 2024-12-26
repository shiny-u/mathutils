//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_CONSTANTS_H
#define MATHUTILS_CONSTANTS_H

#include <cmath>

#ifndef M_PI
    #define M_PI       3.14159265358979323846
#endif

#define MU_JJ std::complex<double>(0., 1.)


#define MU_SQRT2 M_SQRT2
#define MU_SQRT2_2 (MU_SQRT2/2.)

// TODO: voir pour l'integration du framework de gestion d'unites : github.com/nholthaus/units
// Dans ce cas, les constantes relatives aux unites devront etre supprimees


// Time related constants
#define MU_ONE_MINUTE 60.                        ///> NUMBER OF SECONDS IN ONE MINUTE
#define MU_ONE_HOUR (MU_ONE_MINUTE*60.)           ///> NUMBER OF SECONDS IN ONE HOUR
#define MU_ONE_DAY (MU_ONE_HOUR*24.)              ///> NUMBER OF SECONDS UN ONE DAY

// Distance related constants
#define MU_ONE_MILE 1852.                        ///> NUMBER OF METER IN ONE NAUTICAL MILE
#define MU_ONE_KM 1000.                          ///> ONE KILOMETER

// Velocity related constants
#define MU_KNOT (MU_ONE_MILE/MU_ONE_HOUR)          ///> Conversion coeff knot -> m/s
#define MU_KMH (MU_ONE_KM/MU_ONE_HOUR)             ///> ONE KILOMETER BY HOUR

// Angles related constants
#define MU_PI M_PI
#define MU_PI_180 (MU_PI/180.)                    ///> Conversion DEG->RAD
#define DEG2RAD MU_PI_180
#define MU_180_PI (180./MU_PI)                    ///< Conversion RAD->DEG
#define RAD2DEG MU_180_PI
//#define MU_RAD MU_PI_180                          ///> Conversion DEG->RAD
//#define MU_DEG (180./MU_PI)
#define MU_2PI (2.*MU_PI)                         ///> 2*PI
#define MU_PI_2 (MU_PI*0.5)                       ///> PI/2.
#define MU_2_PI (1./MU_PI_2)                      ///> 2./PI
#define MU_4PI (4.*MU_PI)                         ///> 4*PI
#define MU_PI_4 (MU_PI*0.25)                      ///> PI/4
#define MU_1_4PI (1/MU_4PI)                       ///> 1 / (4*PI)

#define MU_RPM2RADS (MU_2PI/60.)                  ///< 1 round per minute is 2*pi rad per 60 seconds

#endif //MATHUTILS_CONSTANTS_H
