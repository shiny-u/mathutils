//
// Created by frongere on 22/12/17.
//

#include "MathUtils/MathUtils.h"
//#include "matplotlibcpp.h"

#include <gtest/gtest.h>


using namespace mathutils;


template <class Scalar>
std::vector<Scalar> MakeData(const std::vector<Scalar> &timeVector, double freq) {

    auto len = timeVector.size();

    std::vector<Scalar> vector;
    vector.reserve(len);

    for (auto t: timeVector) {
        vector.push_back(
                  1.
                + 1 * sin(MU_2PI * freq * t)
                + 2 * sin(MU_2PI * 2*freq * t)
                + 3 * cos(MU_2PI * 5.2*freq * t)
        );
    }

    return vector;

}


TEST(Fft, OldMain) {
//    auto y = Heaviside<double>(-20, 20, 0, 0.01);

//    matplotlibcpp::plot(y);
//    matplotlibcpp::show();


    // Building data
//    unsigned int nfft = Pow2(15);
    unsigned int nfft = Pow2(9);
    std::cout << "NFFT = " << nfft << std::endl;

    auto time = linspace<double>(0., 10., nfft);
    double fs = 1 / (time[1]-time[0]);
    std::cout << "Fs = " << fs << " Hz" << std::endl;

    auto timeVector = MakeData<double>(time, 1);


//    matplotlibcpp::plot(time, timeVector);
//    matplotlibcpp::show();


    // Building the fft object
    FFT<double> fft;
    fft.ScalingON();
//    fft.ScalingOFF();
    fft.HalfSpectrumON();
//    fft.HalfSpectrumOFF();

    // Computing FFT
    std::vector<std::complex<double>> freqVect;
    std::vector<double> frequencies;
    fft.fft(freqVect, frequencies, timeVector, fs, HERTZ);

    // Ploting amplitude and phase spectrums
//    matplotlibcpp::subplot(2, 1, 1);
//    matplotlibcpp::plot(frequencies, Amplitude(freqVect));
//    matplotlibcpp::grid(true);
//    matplotlibcpp::subplot(2, 1, 2);
//    matplotlibcpp::plot(frequencies, Phase(freqVect, DEG));
//    matplotlibcpp::grid(true);
//    matplotlibcpp::show();


    // Computing IFFT from spectrum
    std::vector<double> timeVectRebuild;
    fft.ifft(timeVectRebuild, freqVect);

    EXPECT_TRUE(IsClose(timeVector, timeVectRebuild));
}
