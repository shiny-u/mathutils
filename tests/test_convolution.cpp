//
// Created by frongere on 29/12/17.
//

#include "MathUtils/MathUtilsBoost.h"
//#include "matplotlibcpp.h"

//#include "boost/circular_buffer.hpp"

#include <gtest/gtest.h>


using namespace mathutils;


// TODO delete this unused class
template <class Real>
void ZeroPad(std::vector<Real>& vector, unsigned int n) {
    assert(n >= vector.size());

    vector.reserve(n);
    for (unsigned int i=0; i<vector.size()-n; i++) {
        vector.push_back(0.);
    }
}


TEST(Convolution, OldMain) {
    // TODO Add assertions

    unsigned int n = 500;
    unsigned int m = 500;

    // Creating the impulse response signal
    auto tr = linspace<double>(0, 100., n);
    auto dt = tr[1] - tr[0];
    auto kernel = Sinc(tr);

    // Trailing zéros
    for (unsigned int i=200; i<n; i++) {
        kernel[i] = 0.;
    }

//    // Creating the convolution engine
//    Convolution<double> conv(dt, r);

    // Creating the signal into a circular buffer
    auto ts = linspace<double>(0., 100., m);

    // Implementing and populating a circular buffer in reverse order
    boost::circular_buffer<double> signal(m, m, 0.); // Circular buffer initialized with capacity m and m values set to 0 (buffer is full)
    for (unsigned int i=0; i<m; i++) {
        signal.push_front(
                sin(MU_2PI * 0.05 * ts[i]) +
                sin(MU_2PI * 2*0.05 * ts[i]) +
                sin(MU_2PI * 3.5*0.05 * ts[i])
        ); // push_front put the new values in front on the buffer (not push_bask !)
    }

    // Vector construction for visualization
    std::vector<double> sv;
    for (unsigned int i=0; i<signal.size(); i++) {
        sv.push_back(signal[m-i-1]); // L'indice reflete le fait qu'on prend dans l'ordre inverse...
    }

    // Applying the naive convolution
    auto conv = NaiveConvolution(kernel, signal, dt, true);




//    matplotlibcpp::plot(ts, sv);
//    matplotlibcpp::show();
}
