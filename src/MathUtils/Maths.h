//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_MATHS_H
#define MATHUTILS_MATHS_H

namespace mathutils {

    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template <class Scalar>
    inline Scalar Saturate(Scalar val, Scalar valMin, Scalar valMax) {
        double valout = (val > valMax) ?  valMax : val;
        valout  = (valout<valMin) ? valMin : valout;
        return valout;
    }

    template <class Scalar>
    inline Scalar SgnSqrt(Scalar val) {
        return sgn(val) * sqrt(fabs(val));
    }




}  // end namespace mathutils


#endif //MATHUTILS_MATHS_H
