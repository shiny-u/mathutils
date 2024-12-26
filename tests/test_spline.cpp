//
// Created by frongere on 29/12/17.
//

#include "MathUtils/MathUtils.h"
//#include "matplotlibcpp.h"

#include <gtest/gtest.h>


using namespace mathutils;


TEST(Spline, Heaviside) {
#ifdef SKIP_LONG_TESTS
    GTEST_SKIP() << "Skipped because too long.";
#endif
    int n = 1000;
    auto x = linspace<double>(0, 100, n);

    auto y = Heaviside<double>(x, 20);
//    auto y = Door<double>(x, 25, 40, 5.5);
//    auto y = Dirac<double>(x, 10);
//
//    auto y = Sinc(x);

    Spline<double, 1, 3> s(x, y);

    EXPECT_NEAR(s(19.85), -0.06524310, 1e-8);
    EXPECT_NEAR(s(19.9), -0.08633342, 1e-8);
    EXPECT_NEAR(s(20.0), 0.83513843, 1e-8);
    EXPECT_NEAR(s(20.1), 1.04226281, 1e-8);
    EXPECT_NEAR(s(20.15), 0.97220185, 1e-8);
    EXPECT_NEAR(s(21.0), 0.99999968, 1e-8);
    EXPECT_NEAR(s(21.05), 1.00000019, 1e-8);

    // TODO: ajouter methode pour l'evaluation de la spline sur un echantillopn std::vector...
//    int ni = 10000;
//    auto xi = linspace<double>(0, 100, ni);
//    std::vector<double> yi(ni);
//    for (unsigned int i=0; i<ni; i++) {
//        yi[i] = s(xi[i]);
//    }

//    matplotlibcpp::plot(x, y);
//    matplotlibcpp::plot(xi, yi);
//    matplotlibcpp::show();
}
