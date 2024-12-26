//
// Created by pierre-yves on 01/12/2020.
// This test checks the computation of the triple Chebyshev series.
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
class Function3dTest : public Function3d<double> {
public:
    /// This method evaluates the function at the point (x, y, z).
    double Evaluate(const double &x, const double &y, const double &z) const override {
        return 1. / (2. * x + 3. * y + 4. * z);
    }

    /// This method evaluates the x-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_x(const double &x, const double &y, const double &z) const {
        return -2. / pow(2. * x + 3. * y + 4. * z, 2.);
    }

    /// This method evaluates the y-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_y(const double &x, const double &y, const double &z) const {
        return -3. / pow(2. * x + 3. * y + 4. * z, 2.);
    }

    /// This method evaluates the z-derivative of the function at the point (x, y, z).
    double Evaluate_derivative_z(const double &x, const double &y, const double &z) const {
        return -4. / pow(2. * x + 3. * y + 4. * z, 2.);
    }
};


TEST(ChebyshevSeries3d, ClosedSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                              Closed segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function3dTest myFunction3d;
    double xmin = 1.;
    double xmax = 10.;
    double ymin = 1.;
    double ymax = 10.;
    double zmin = 1.;
    double zmax = 10.;
    int order_x = 11;
    int order_y = 13;
    int order_z = 10;
    auto myChebyshevApprox3dClosed = ChebyshevSeries3dClosed<double>(&myFunction3d, xmin, xmax, ymin, ymax, zmin, zmax,
                                                                     order_x, order_y, order_z);
    myChebyshevApprox3dClosed.Compute_aijk();

    // Tests.
    PrintHeader("Closed segments - Function - Zero of a ChebyShev polynomial.");
    int r = 3;
    double tmp = cos(MU_PI_2 * (2. * r + 1.) / (order_x + 1.));
    double x = 0.5 * (xmax - xmin) * tmp + 0.5 * (xmax + xmin);
    int s = 3;
    tmp = cos(MU_PI_2 * (2. * s + 1.) / (order_y + 1.));
    double y = 0.5 * (ymax - ymin) * tmp + 0.5 * (ymax + ymin);
    int t = 3;
    tmp = cos(MU_PI_2 * (2. * t + 1.) / (order_z + 1.));
    double z = 0.5 * (zmax - zmin) * tmp + 0.5 * (zmax + zmin);
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    double num = myChebyshevApprox3dClosed.Evaluate(x, y, z);
    std::cout << "Triple Chebyshev series approximation: " << num << std::endl;
    assert(IsClose(myChebyshevApprox3dClosed.Evaluate(x, y, z), myFunction3d.Evaluate(x, y, z)));

    PrintHeader("Closed segments - Function - Arbitrary point.");
    x = 3.;
    y = 2.;
    z = 5.;
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dClosed.Evaluate(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Closed segments - x derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Closed segments - y derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("Closed segments - z derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_z(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(ChebyshevSeries3d, MixedHalfOpenAndClosedSegments) {
#ifdef SKIP_LONG_TESTS
    GTEST_SKIP() << "Skipped because too long.";
#endif
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                          Mixed half-open and closed segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function3dTest myFunction3d;
    double xmin = 1.;
    double xmax = 30;
    double ymin = 1.;
    double zmin = 1.;
    double order_x = 15;
    double order_y = 30;
    double order_z = 30;
    double x = 10;
    double y = 45;
    double z = 50;
    auto myChebyshevApprox3dYZOpenedXClosed = ChebyshevSeries3dYZOpenedXClosed<double>(&myFunction3d, xmin, xmax, ymin,
                                                                                       zmin, order_x, order_y, order_z);
    myChebyshevApprox3dYZOpenedXClosed.Compute_aijk();

    // Tests.
    PrintHeader("X closed, YZ open - Function - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    double num = myChebyshevApprox3dYZOpenedXClosed.Evaluate(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("X closed, YZ open - x derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dYZOpenedXClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("X closed, YZ open - y derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dYZOpenedXClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("X closed, YZ open - z derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_z(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dYZOpenedXClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(ChebyshevSeries3d, MixedHalfOpenAndClosedSegmentsBis) {
#ifdef SKIP_LONG_TESTS
    GTEST_SKIP() << "Skipped because too long.";
#endif
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                          Mixed half-open and closed segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function3dTest myFunction3d;
    double xmin = 1.;
    double xmax = 30;
    double ymin = 1.;
    double ymax = 30.;
    double zmin = 1.;
    double order_x = 15;
    double order_y = 15;
    double order_z = 30;
    double x = 10;
    double y = 8;
    double z = 50;
    auto myChebyshevApprox3dZOpenedXYClosed = ChebyshevSeries3dZOpenedXYClosed<double>(&myFunction3d, xmin, xmax, ymin,
                                                                                       ymax, zmin, order_x, order_y,
                                                                                       order_z);
    myChebyshevApprox3dZOpenedXYClosed.Compute_aijk();

    // Tests.
    PrintHeader("XY closed, Z open - Function - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    double num = myChebyshevApprox3dZOpenedXYClosed.Evaluate(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("XY closed, Z open - x derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dZOpenedXYClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("XY closed, Z open - y derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dZOpenedXYClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));

    PrintHeader("XY closed, Z open - z derivative - Arbitrary point.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_z(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    num = myChebyshevApprox3dZOpenedXYClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << num << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((num - ana) / ana) << std::endl;
    assert(IsClose(num, ana, 10e-3));
}
