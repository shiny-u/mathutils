//
// Created by pierre-yves on 04/03/2021.
// This test checks the transformation from double Chebyshev series to double power series.
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


void PrintInfo(std::string info) {
  std::cout << info << ":" << std::endl;
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
        return -2. / pow(2. * x + 3. * y, 2.);
    }

    /// This method evaluates the y-derivative of the function at the point (x, y).
    double Evaluate_derivative_y(const double &x, const double &y) const {
        return -3. / pow(2. * x + 3. * y, 2.);
    }
};


TEST(PowerSeries2d, ClosedSegments) {
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
    auto myChebyshevSeries2dClosed = ChebyshevSeries2dClosed<double>(&myFunction2d, xmin, xmax, ymin, ymax, order_x,
                                                                     order_y);
    myChebyshevSeries2dClosed.Compute_aij();

    // From Chebyshev series to power series.
    auto bij = myChebyshevSeries2dClosed.Compute_bij();
    auto myPowerSeries2dClosed = PowerSeries2dClosed<double>(bij, xmin, xmax, ymin, ymax);

    PrintHeader("Closed segments - Function.");
    double x = 3.;
    double y = 2.;
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    double ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries2dClosed.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries2dClosed.Evaluate(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Closed segments - x derivative.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_x(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries2dClosed.Evaluate_derivative_x(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries2dClosed.Evaluate_derivative_x(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Closed segments - y derivative.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_y(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries2dClosed.Evaluate_derivative_y(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries2dClosed.Evaluate_derivative_y(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(PowerSeries2d, HalfOpenSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                            Half-open segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function2dTest myFunction2d;
    double xmin = 1.;
    double ymin = 1.;
    double order_x = 18;
    double order_y = 18;
    double x = 15;
    double y = 16;
    auto myChebyshevSeries2dOpened = ChebyshevSeries2dOpened<double>(&myFunction2d, xmin, ymin, order_x, order_y);
    myChebyshevSeries2dOpened.Compute_aij();

    // From Chebyshev series to power series.
    auto bij = myChebyshevSeries2dOpened.Compute_bij();
    auto myPowerSeries2dOpened = PowerSeries2dOpened<double>(bij, xmin, ymin);

    // Tests.
    PrintHeader("Half-open segments - Function.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    double ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries2dOpened.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries2dOpened.Evaluate(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Half-open segments - x derivative.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_x(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries2dOpened.Evaluate_derivative_x(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries2dOpened.Evaluate_derivative_x(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Half-open segments - y derivative.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_y(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries2dOpened.Evaluate_derivative_y(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries2dOpened.Evaluate_derivative_y(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(PowerSeries2d, MixedClosedHalfOpenSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                            Mixed closed - half-open segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function2dTest myFunction2d;
    double xmin = 1.;
    double xmax = 10.;
    double ymin = 1.;
    double order_x = 10;
    double order_y = 18;
    double x = 3;
    double y = 4;
    auto myChebyshevSeries2dMixed = ChebyshevSeries2dMixed<double>(&myFunction2d, xmin, xmax, ymin, order_x, order_y);
    myChebyshevSeries2dMixed.Compute_aij();

    // From Chebyshev series to power series.
    auto bij = myChebyshevSeries2dMixed.Compute_bij();
    auto myPowerSeries2dMixed = PowerSeries2dMixed<double>(bij, xmin, xmax, ymin);

    // Tests.
    PrintHeader("Mixed segments - Function.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    double ana = myFunction2d.Evaluate(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries2dMixed.Evaluate(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries2dMixed.Evaluate(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Mixed segments - x derivative.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_x(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries2dMixed.Evaluate_derivative_x(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries2dMixed.Evaluate_derivative_x(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Mixed segments - y derivative.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction2d.Evaluate_derivative_y(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries2dMixed.Evaluate_derivative_y(x, y);
    std::cout << "Double Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries2dMixed.Evaluate_derivative_y(x, y);
    std::cout << "Double Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));
}
