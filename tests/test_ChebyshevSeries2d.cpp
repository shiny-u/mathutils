//
// Created by pierre-yves on 01/12/2020.
// This test checks the computation of the double Chebyshev series.
//

#include <iostream>
#include "MathUtils/MathUtilsBoost.h"

#include <gtest/gtest.h>


using namespace mathutils;


void PrintHeader(std::string title) {
  std::cout << "\n=====================================================================" << std::endl;
  std::cout << "    " << title << std::endl;
  std::cout << "=====================================================================" << std::endl;
}


// Definition of a function to approximate by a double Chebyshev series.
class Function2dTest : public Function2d<double> {
public:
    /// This method evaluates the function at the point (x, y).
    double Evaluate(const double &x, const double &y) const override {
        return 1. / (2. * x + 3. * y);
    }
    /// This method evaluates the x-derivative of the function at the point (x, y).
    double Evaluate_derivative_x(const double &x, const double &y) const {
        return - 2. / pow(2. * x + 3. * y, 2.);
    }
    /// This method evaluates the y-derivative of the function at the point (x, y).
    double Evaluate_derivative_y(const double &x, const double &y) const {
        return - 3. / pow(2. * x + 3. * y, 2.);
    }
};


void PrintInfo(std::string info) {
  std::cout << info << ":" << std::endl;
}


TEST(ChebyshevSeries2d, ClosedSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                              Closed segments.
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function2dTest myFunction2d;
    double xmin = 1.;
    double xmax = 10.;
    double ymin = 1.;
    double ymax = 10.;
    int order_x = 11;
    int order_y = 13;
    auto myChebyshevApprox2dClosed = ChebyshevSeries2dClosed<double>(&myFunction2d, xmin, xmax, ymin, ymax, order_x,
                                                                     order_y);
    myChebyshevApprox2dClosed.Compute_aij();

    // Tests.
    PrintHeader("Closed segments - Function - Zero of a ChebyShev polynomial.");
    int r = 3;
    double tmp = cos(MU_PI_2 * (2. * r + 1.) / (order_x + 1.));
    double x = 0.5 * (xmax - xmin) * tmp + 0.5 * (xmax + xmin);
    int s = 3;
    tmp = cos(MU_PI_2 * (2. * s + 1.) / (order_y + 1.));
    double y = 0.5 * (ymax - ymin) * tmp + 0.5 * (ymax + ymin);
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    double ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    double num = myChebyshevApprox2dClosed.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    assert(IsClose(num, ana));

    PrintHeader("Closed segments - Function - Arbitrary point.");
    x = 3.;
    y = 2.;
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dClosed.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Closed segments - x derivative - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_x(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dClosed.Evaluate_derivative_x(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Closed segments - y derivative - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_y(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dClosed.Evaluate_derivative_y(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(ChebyshevSeries2d, HalfOpenSegments) {
#ifdef SKIP_LONG_TESTS
    GTEST_SKIP() << "Skipped because too long.";
#endif
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                            Half-open segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function2dTest myFunction2d;
    double xmin = 1.;
    double ymin = 1.;
    double order_x = 50;
    double order_y = 50;
    double x = 130;
    double y = 150;
    auto myChebyshevApprox2dOpened = ChebyshevSeries2dOpened<double>(&myFunction2d, xmin, ymin, order_x, order_y);
    myChebyshevApprox2dOpened.Compute_aij();

    // Tests.
    PrintHeader("Half-open segments - Function - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    double ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    double num = myChebyshevApprox2dOpened.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Half-open segments - x derivative - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_x(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dOpened.Evaluate_derivative_x(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Half-open segments - y derivative - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_y(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dOpened.Evaluate_derivative_y(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(ChebyshevSeries2d, MixedClosedHalfOpenSegments) {
#ifdef SKIP_LONG_TESTS
    GTEST_SKIP() << "Skipped because too long.";
#endif
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                            Mixed closed - half-open segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function2dTest myFunction2d;
    double xmin = 1.;
    double xmax = 10.;
    double ymin = 1.;
    double order_x = 50;
    double order_y = 50;
    double x = 3;
    double y = 150;
    auto myChebyshevApprox2dMixed = ChebyshevSeries2dMixed<double>(&myFunction2d, xmin, xmax, ymin, order_x, order_y);
    myChebyshevApprox2dMixed.Compute_aij();

    // Tests.
    PrintHeader("Mixed segments - Function - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    double ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    double num = myChebyshevApprox2dMixed.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Mixed segments - x derivative - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_x(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dMixed.Evaluate_derivative_x(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Mixed segments - y derivative - Arbitrary point.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_y(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox2dMixed.Evaluate_derivative_y(x, y);
    std::cout << "Double Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));
}
