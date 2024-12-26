//
// Created by pierre-yves on 01/12/2020.
// This test suite checks the computation of the single Chebyshev series.
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
        return - 0.5 / (x * x);
    }
};


TEST(ChebyshevSeries1d, ClosedSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                              Closed segments.
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function1dTest myFunction1d;
    double xmin = 1.;
    double xmax = 10.;
    int order_x = 11;
    auto myChebyshevApprox1dClosed = ChebyshevSeries1dClosed<double>(&myFunction1d, xmin, xmax, order_x);
    myChebyshevApprox1dClosed.Compute_ai();

    // Tests.
    {
        PrintHeader("Closed segments - Function - Zero of a ChebyShev polynomial.");
        int r = 3;
        double tmp = cos(MU_PI_2 * (2. * r + 1.) / (order_x + 1.));
        double x = 0.5 * (xmax - xmin) * tmp + 0.5 * (xmax + xmin);
        std::cout << "Point: x = " << x << std::endl;
        double ana = myFunction1d.Evaluate(x);
        std::cout << "Analytical: " << ana << std::endl;
        double num = myChebyshevApprox1dClosed.Evaluate(x);
        std::cout << "Single Chebyshev series approximation: " << num << std::endl;
        EXPECT_NEAR(num, ana, 1e-8);
    }

    {
        PrintHeader("Closed segments - Function - Arbitrary point.");
        double x = 3.;
        std::cout << "Point: x = " << x << std::endl;
        double ana = myFunction1d.Evaluate(x);
        std::cout << "Analytical: " << ana << std::endl;
        double num = myChebyshevApprox1dClosed.Evaluate(x);
        std::cout << "Single Chebyshev series approximation: " << num << std::endl;
        std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
        EXPECT_NEAR(num, ana, 1e-2);
    }

    {
        PrintHeader("Closed segments - x derivative - Arbitrary point.");
        double x = 3.;
        std::cout << "Point: x = " << x << std::endl;
        double ana = myFunction1d.Evaluate_derivative_x(x);
        std::cout << "Analytical: " << ana << std::endl;
        double num = myChebyshevApprox1dClosed.Evaluate_derivative_x(x);
        std::cout << "Single Chebyshev series approximation: " << num << std::endl;
        std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
        EXPECT_NEAR(num, ana, 1e-2);
    }

    std::cout << "\n" << std::endl;
}


TEST(ChebyshevSeries1d, HalfOpenSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                            Half-open segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function1dTest myFunction1d;
    double xmin = 1.;
    int order_x = 50;
    double x = 150;
    auto myChebyshevApprox1dOpened = ChebyshevSeries1dOpened<double>(&myFunction1d, xmin, order_x);
    myChebyshevApprox1dOpened.Compute_ai();

    // Tests.
    {
        PrintHeader("Half-open segments - Function - Arbitrary point.");
        std::cout << "Point: x = " << x << std::endl;
        double ana = myFunction1d.Evaluate(x);
        std::cout << "Analytical: " << ana << std::endl;
        double num = myChebyshevApprox1dOpened.Evaluate(x);
        std::cout << "Single Chebyshev series approximation: " << num << std::endl;
        std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
        EXPECT_NEAR(num, ana, 1e-2);
    }

    {
        PrintHeader("Half-open segments - x derivative - Arbitrary point.");
        std::cout << "Point: x = " << x << std::endl;
        double ana = myFunction1d.Evaluate_derivative_x(x);
        std::cout << "Analytical: " << ana << std::endl;
        double num = myChebyshevApprox1dOpened.Evaluate_derivative_x(x);
        std::cout << "Single Chebyshev series approximation: " << num << std::endl;
        std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
        EXPECT_NEAR(num, ana, 1e-2);
    }
}
