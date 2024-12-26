//
// Created by pierre-yves on 04/03/2021.
// This test checks the transformation from single Chebyshev series to single power series.
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


// Definition of a function to approximate by a single Chebyshev series.
class Function1dTest : public Function1d<double> {
public:
    /// This method evaluates the function at the point x.
    double Evaluate(const double &x) const override {
        return 1. / (2. * x);
    }

    /// This method evaluates the x-derivative of the function at the point x.
    double Evaluate_derivative_x(const double &x) const {
        return -0.5 / (x * x);
    }

};


TEST(PowerSeries1d, ClosedSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                              Closed segments.
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ai coefficients.
    Function1dTest myFunction1d;
    double xmin = 1.;
    double xmax = 10.;
    int order_x = 11;
    auto myChebyshevSeries1dClosed = ChebyshevSeries1dClosed<double>(&myFunction1d, xmin, xmax, order_x);
    myChebyshevSeries1dClosed.Compute_ai();

    // From Chebyshev series to power series.
    auto bi = myChebyshevSeries1dClosed.Compute_bi();
    auto myPowerSeries1dClosed = PowerSeries1dClosed<double>(bi, xmin, xmax);

    PrintHeader("Closed segments - Function.");
    double x = 3.;
    std::cout << "Point: x = " << x << std::endl;
    double ana = myFunction1d.Evaluate(x);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries1dClosed.Evaluate(x);
    std::cout << "Single Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries1dClosed.Evaluate(x);
    std::cout << "Single Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Closed segments - x derivative.");
    std::cout << "Point: x = " << x << std::endl;
    ana = myFunction1d.Evaluate_derivative_x(x);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries1dClosed.Evaluate_derivative_x(x);
    std::cout << "Single Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries1dClosed.Evaluate_derivative_x(x);
    std::cout << "Single Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    std::cout << "\n" << std::endl;
}

TEST(PowerSeries1d, HalfOpenSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                            Half-open segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function1dTest myFunction1d;
    double xmin = 1.;
    double order_x = 18;
    double x = 15;
    auto myChebyshevSeries1dOpened = ChebyshevSeries1dOpened<double>(&myFunction1d, xmin, order_x);
    myChebyshevSeries1dOpened.Compute_ai();

    // From Chebyshev series to power series.
    auto bi = myChebyshevSeries1dOpened.Compute_bi();
    auto myPowerSeries1dOpened = PowerSeries1dOpened<double>(bi, xmin);

    // Tests.
    PrintHeader("Half-open segments - Function.");
    std::cout << "Point: x = " << x << std::endl;
    double ana = myFunction1d.Evaluate(x);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries1dOpened.Evaluate(x);
    std::cout << "Single Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries1dOpened.Evaluate(x);
    std::cout << "Single Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Half-open segments - x derivative.");
    std::cout << "Point: x = " << x << std::endl;
    ana = myFunction1d.Evaluate_derivative_x(x);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries1dOpened.Evaluate_derivative_x(x);
    std::cout << "Single Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries1dOpened.Evaluate_derivative_x(x);
    std::cout << "Single Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));
}
