//
// Created by frongere on 22/12/17.
//

#ifndef MATHUTILS_FFT_H
#define MATHUTILS_FFT_H

//#include "kiss_fft.h"

#include <unsupported/Eigen/FFT>
#include "Unit.h"

namespace mathutils {

    template <class Scalar>
    class FFT : protected Eigen::FFT<Scalar> {

    public:

        // TODO: voir pour un zero padding auto ?
        void fft(std::vector<std::complex<Scalar>>& freqVector, const std::vector<Scalar>& vector) {
            // TODO: voir le comportement si on ne donne pas vector en puissance de 2...
            this->fwd(freqVector, vector);
        }

        void fft(std::vector<std::complex<Scalar>>& freqVector, std::vector<Scalar>& frequencies,
                 const std::vector<Scalar>& vector, const Scalar sampleFrequency, FREQUENCY_UNIT frequencyUnit=HERTZ) {
            // TODO: voir le comportement si on ne donne pas vector en puissance de 2...
            // TODO: fournir plutot le vecteur temps et e deduire directement la frequence d'echantilonnage
            this->fwd(freqVector, vector);
            frequencies = GetFrequencyVector((unsigned int)vector.size(), sampleFrequency, frequencyUnit);
        }

        void ifft(std::vector<Scalar>& vector, const std::vector<std::complex<Scalar>>& freqVector) {
            this->inv(vector, freqVector);
        }

        void HalfSpectrumON() {
            this->SetFlag(Eigen::FFT<Scalar>::HalfSpectrum);
//            this->HalfSpectrum = true;
        }

        void HalfSpectrumOFF() {
            this->ClearFlag(Eigen::FFT<Scalar>::HalfSpectrum);
//            this->HalfSpectrum = false;
        }

        void ScalingON() {
            this->ClearFlag(Eigen::FFT<Scalar>::Unscaled);
//            this->Unscaled = false;
        }

        void ScalingOFF() {
            this->SetFlag(Eigen::FFT<Scalar>::Unscaled);
//            this->Unscaled = true;
        }

        std::vector<Scalar> GetFrequencyVector(const unsigned int nfft, const Scalar sampleFrequency,
                                               FREQUENCY_UNIT frequencyUnit=HERTZ) {

            auto df = sampleFrequency / (Scalar)nfft;

            unsigned int n = nfft;
            if (this->HasFlag(Eigen::FFT<Scalar>::HalfSpectrum)) {
                n = n/2 + 1;
            }

            std::vector<Scalar> out(n);

            Scalar f = 0.;
            for (unsigned int i=0; i<n; i++) {
                out[i] = convert_frequency(f, HERTZ, frequencyUnit);
                f += df;
            }

            return out;

        }



    };


    template <class Real>
    std::vector<Real> ConvolveNaive(const std::vector<Real>& f, const std::vector<Real>& g) {
        auto nf = f.size();
        auto ng = g.size();

        auto nc = nf + g -1;




    }



// FIXME: voir pourquoi NextPow2 est rejete des lors qu'on a FFT.h dans FRyDoM !!! (definition multiple...)

//    unsigned int NextPow2(const unsigned int n, const bool returnValue=true) {
//        unsigned int p = 0;
//        unsigned int val;
//        while (true) {
//            val = (unsigned int)pow(2, p);
//            if (val >= n) {
//                if (returnValue) {
//                    return val;
//                } else {
//                    return p;
//                }
//            }
//            p++;
//        }
//    }

    inline unsigned int Pow2(unsigned int p) {
        return 1 << p;
    }



}  // end namespace mathutils


#endif //MATHUTILS_FFT_H
