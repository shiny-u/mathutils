//
// Created by frongere on 29/12/17.
//

#ifndef MATHUTILS_CONVOLUTION_H
#define MATHUTILS_CONVOLUTION_H

#include "boost/circular_buffer.hpp"
//#include "matplotlibcpp.h"

namespace mathutils {


    template <class Real>
    class Convolution {

    private:
        Real m_timeStep;
        std::vector<Real> m_kernel;

    public:
        Convolution() = default;

        explicit Convolution(const Real dt, const std::vector<Real>& kernel)
                : m_timeStep(dt), m_kernel(kernel) {}

        void SetKernel(const Real& timeStep, const std::vector<Real>& kernel) {
            m_timeStep = timeStep;
            m_kernel = kernel;
        }



    };



    // TODO: creer une classe dans laquelle on peut renseigner le noyau de convolution (reponse impulsionnelle) et fournir un signal
    // TODO: attention, on veut avoir un kernel et un signal echantillonnes de la meme maniere !!!
    template <class Real>
    Real NaiveConvolution(const std::vector<Real>& kernel, const boost::circular_buffer<Real>& signal, Real dt, bool reverse=true) {

        unsigned long N = kernel.size();

        assert(signal.size() >= N);

        // Commputing the integrand
        std::vector<Real> integrand(N);
        if (reverse) {
            for (unsigned long i = 0; i < N; i++) { // TODO: gerer le cas ou reverse est false (il faut renverser le signal ici ou creer une boucle differente)
                integrand[i] = kernel[i] * signal[i];
            }
        } else {
            for (unsigned long i = 0; i < N; i++) { // TODO: gerer le cas ou reverse est false (il faut renverser le signal ici ou creer une boucle differente)
                integrand[i] = kernel[i] * signal[N-i-1];
            }
        }

        // A retirer, test
        std::vector<double> sv;
        for (unsigned int i=0; i<signal.size(); i++) {
            sv.push_back(signal[i]); // L'indice reflete le fait qu'on prend dans l'ordre inverse...
        }


//        matplotlibcpp::plot(kernel);
//        matplotlibcpp::plot(sv);
//        matplotlibcpp::plot(integrand);
//        matplotlibcpp::show();


        // Computing the integral (with the trapezoidal rule for instance)
        return Trapz(integrand, dt);

    }


}  // end namespace mathutils

#endif //MATHUTILS_CONVOLUTION_H
