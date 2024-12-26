//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_VECTORGENERATION_H
#define MATHUTILS_VECTORGENERATION_H

#include <cassert>
#include <cmath>
#include <vector>


namespace mathutils {

    template<class Real=double>
    std::vector<Real> linspace(Real start, Real stop, unsigned int num) {
        assert(num > 1);
        assert(start <= stop);

        unsigned int numm1 = num-1;

        Real step = (stop-start) / (numm1);

        std::vector<Real> out(num);
        for (unsigned long i=0; i<numm1; ++i) {
            out[i] = start + i*step;
        }
        out[numm1] = stop;
        return out;
    }

    template <class Real=double>
    std::vector<Real> logspace(Real start, Real stop, unsigned int num, Real base=10.) {
        assert(num > 1);
        assert(start <= stop);

        auto out = linspace(start, stop, num);
        for (unsigned int i = 0; i < out.size(); i++) {
            out.at(i) = pow(base, out.at(i));
        }
        return out;
    }

    template<class Real=double>
    std::vector<Real> arange(Real start, Real stop, Real step=1) {
        assert (start < stop);
        assert (step > 0.);

        auto num = floor((stop-start)/step);
        auto final_stop = num * step;

        return linspace(start, final_stop, (unsigned int)num);
    }


}


#endif //MATHUTILS_VECTORGENERATION_H
