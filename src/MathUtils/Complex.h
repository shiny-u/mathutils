//
// Created by frongere on 28/12/17.
//

#ifndef MATHUTILS_COMPLEX_H
#define MATHUTILS_COMPLEX_H

namespace mathutils {

    template <class Real_t>
    std::vector<Real_t> Real(const std::vector<std::complex<Real_t>>& vector) {
        unsigned long n = vector.size();
        std::vector<Real_t> out(n);
        for (unsigned long i=0; i<n; i++) {
            out[i] = vector[i].real();
        }
        return out;
    }

    template <class Real_t>
    std::vector<Real_t> Imag(const std::vector<std::complex<Real_t>>& vector) {
        unsigned long n = vector.size();
        std::vector<Real_t> out(n);
        for (unsigned long i=0; i<n; i++) {
            out[i] = vector[i].imag();
        }
        return out;
    }

    template <class Real_t>
    std::vector<Real_t> Amplitude(const std::vector<std::complex<Real_t>>& vector) {
        unsigned long n = vector.size();
        std::vector<Real_t> out(n);
        Real_t x, y;
        for (unsigned long i=0; i<n; i++) {
            x = vector[i].real();
            y = vector[i].imag();
            out[i] = sqrt(x*x + y*y);
        }
        return out;
    }

    template <class Real_t>
    std::vector<Real_t> Phase(const std::vector<std::complex<Real_t>>& vector, ANGLE_UNIT unit=RAD) {

        Real_t c = 1.;
        if (unit == DEG) {
            c *= MU_180_PI;
        }

        unsigned long n = vector.size();
        std::vector<Real_t> out(n);
        Real_t x, y;
        for (unsigned long i=0; i<n; i++) {
            x = vector[i].real();
            y = vector[i].imag();
            out[i] = atan2(y, x) * c;
//            out[i] = atan(y/x) * c;
        }
        return out;
    }




}


#endif //MATHUTILS_COMPLEX_H
