//
// Created by pierre-yves on 05/03/2021.
// This test checks the transformation from triple Chebyshev series to triple power series.
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
        return 1. + 1. / (2. * x + 3. * y + 4. * z);
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

    /// This method evaluates the function at the point (x, y) for z tending to +infinity.
    double Evaluate_zinf(const double &x, const double &y) const {
        return 1.;
    }

    /// This method evaluates the x-derivative of the function at the point (x, y) for z tending to +infinity.
    double Evaluate_derivative_x_zinf(const double &x, const double &y) const {
        return 0.;
    }

    /// This method evaluates the y-derivative of the function at the point (x, y) for z tending to +infinity.
    double Evaluate_derivative_y_zinf(const double &x, const double &y) const {
        return 0.;
    }

};


TEST(PowerSeries3d, ClosedSegments) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                              Closed segments
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function3dTest myFunction3d;
    double xmin = 1.;
    double xmax = 6.;
    double ymin = 1.;
    double ymax = 6.;
    double zmin = 1.;
    double zmax = 6.;
    int order_x = 11;
    int order_y = 13;
    int order_z = 10;
    auto myChebyshevSeries3dClosed = ChebyshevSeries3dClosed<double>(&myFunction3d, xmin, xmax, ymin, ymax, zmin, zmax,
                                                                     order_x, order_y, order_z);
    myChebyshevSeries3dClosed.Compute_aijk();

    // From Chebyshev series to power series.
    auto bijk = myChebyshevSeries3dClosed.Compute_bijk();
    auto myPowerSeries3dClosed = PowerSeries3dClosed<double>(bijk, xmin, xmax, ymin, ymax, zmin, zmax);

    PrintHeader("Closed segments - Function.");
    double x = 3.;
    double y = 2.;
    double z = 5.;
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries3dClosed.Evaluate(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries3dClosed.Evaluate(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Closed segments - x derivative.");
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Closed segments - y derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("Closed segments - z derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_z(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));
    std::cout << "\n" << std::endl;
}


TEST(PowerSeries3d, ClosedSegmentsPredefinedZ) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                          Closed segments - Predefined z
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function3dTest myFunction3d;
    double xmin = 1.;
    double xmax = 6.;
    double ymin = 1.;
    double ymax = 6.;
    double zmin = 1.;
    double zmax = 6.;
    int order_x = 11;
    int order_y = 13;
    int order_z = 10;
    auto myChebyshevSeries3dClosed = ChebyshevSeries3dClosed<double>(&myFunction3d, xmin, xmax, ymin, ymax, zmin, zmax,
                                                                     order_x, order_y, order_z);
    myChebyshevSeries3dClosed.Compute_aijk();

    // From Chebyshev series to power series.
    auto bijk = myChebyshevSeries3dClosed.Compute_bijk();
    auto myPowerSeries3dClosed = PowerSeries3dClosed<double>(bijk, xmin, xmax, ymin, ymax, zmin, zmax);

    PrintHeader("Closed segments - Function - Predefined z.");
    double x = 3.;
    double y = 2.;
    double z = 5.;
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    clock_t t_3d = clock();
    double power = myPowerSeries3dClosed.Evaluate(x, y, z);
    t_3d = clock() - t_3d;
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    auto cij = myPowerSeries3dClosed.Compute_cij(z);
    clock_t t_2d = clock();
    double power_predefined_z = myPowerSeries3dClosed.Evaluate_z_predefined(x, y, cij);
    t_2d = clock() - t_2d;
    std::cout << "Triple Power series approximation with predefined z: " << power_predefined_z << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power_predefined_z - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, power_predefined_z, 10e-3));
    std::cout << "CPU time triple power series (s) : " << ((float) t_3d / CLOCKS_PER_SEC) << std::endl;
    std::cout << "CPU time double power series with predefined z (s) : " << ((float) t_2d / CLOCKS_PER_SEC)
              << std::endl;
    std::cout << "\n" << std::endl;

    PrintHeader("Closed segments - x derivative - Predefined z.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    t_3d = clock();
    power = myPowerSeries3dClosed.Evaluate_derivative_x(x, y, z);
    t_3d = clock() - t_3d;
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    t_2d = clock();
    power_predefined_z = myPowerSeries3dClosed.Evaluate_derivative_x_z_predefined(x, y, cij);
    t_2d = clock() - t_2d;
    std::cout << "Triple Power series approximation with predefined z: " << power_predefined_z << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power_predefined_z - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, power_predefined_z, 10e-3));
    std::cout << "CPU time triple power series (s) : " << ((float) t_3d / CLOCKS_PER_SEC) << std::endl;
    std::cout << "CPU time double power series with predefined z (s) : " << ((float) t_2d / CLOCKS_PER_SEC)
              << std::endl;
    std::cout << "\n" << std::endl;

    PrintHeader("Closed segments - y derivative - Predefined z.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    t_3d = clock();
    power = myPowerSeries3dClosed.Evaluate_derivative_y(x, y, z);
    t_3d = clock() - t_3d;
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    t_2d = clock();
    power_predefined_z = myPowerSeries3dClosed.Evaluate_derivative_y_z_predefined(x, y, cij);
    t_2d = clock() - t_2d;
    std::cout << "Triple Power series approximation with predefined z: " << power_predefined_z << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power_predefined_z - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, power_predefined_z, 10e-3));
    std::cout << "CPU time triple power series (s) : " << ((float) t_3d / CLOCKS_PER_SEC) << std::endl;
    std::cout << "CPU time double power series with predefined z (s) : " << ((float) t_2d / CLOCKS_PER_SEC)
              << std::endl;
    std::cout << "\n" << std::endl;
}


TEST(PowerSeries3d, ClosedSegmentsFactorizationPredefinedZ) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                          Closed segments - Factorization - Predefined z
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // aij coefficients.
    Function3dTest myFunction3d;
    double xmin = 1.;
    double xmax = 6.;
    double ymin = 1.;
    double ymax = 6.;
    double zmin = 1.;
    double zmax = 6.;
    int order_x = 11;
    int order_y = 13;
    int order_z = 10;
    auto myChebyshevSeries3dClosed = ChebyshevSeries3dClosed<double>(&myFunction3d, xmin, xmax, ymin, ymax, zmin, zmax,
                                                                     order_x, order_y, order_z);
    myChebyshevSeries3dClosed.Compute_aijk();

    // From Chebyshev series to power series.
    auto bijk = myChebyshevSeries3dClosed.Compute_bijk();
    auto myPowerSeries3dClosed = PowerSeries3dClosed<double>(bijk, xmin, xmax, ymin, ymax, zmin, zmax);

    PrintHeader("Closed segments - Factorization Predefined z.");
    double x = 3.;
    double y = 2.;
    double z = 5.;
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    auto cij = myPowerSeries3dClosed.Compute_cij(z);
    clock_t t_no_factorization = clock();
    double f_no_factorization, dfdx_no_factorization, dfdy_no_factorization;
    myPowerSeries3dClosed.Evaluate_f_dfdx_dfdy(x, y, z, f_no_factorization, dfdx_no_factorization,
                                               dfdy_no_factorization);
    t_no_factorization = clock() - t_no_factorization;
    std::cout << "f, dfdx, dfdy without factorization = " << " " << f_no_factorization << " " << dfdx_no_factorization
              << " " << dfdy_no_factorization << std::endl;
    cij = myPowerSeries3dClosed.Compute_cij(z);
    clock_t t_with_factorization = clock();
    double f_with_factorization, dfdx_with_factorization, dfdy_with_factorization;
    myPowerSeries3dClosed.Evaluate_f_dfdx_dfdy_z_predefined(x, y, cij, f_with_factorization, dfdx_with_factorization,
                                                            dfdy_with_factorization);
    t_with_factorization = clock() - t_with_factorization;
    std::cout << "f, dfdx, dfdy with factorization = " << " " << f_with_factorization << " " << dfdx_with_factorization
              << " " << dfdy_with_factorization << std::endl;
    EXPECT_TRUE(IsClose(f_no_factorization, f_with_factorization, 10e-3));
    EXPECT_TRUE(IsClose(dfdx_no_factorization, dfdx_with_factorization, 10e-3));
    EXPECT_TRUE(IsClose(dfdy_no_factorization, dfdy_with_factorization, 10e-3));
    std::cout << "CPU time without factorization (s) : " << ((float) t_no_factorization / CLOCKS_PER_SEC) << std::endl;
    std::cout << "CPU time with factorization (s) : " << ((float) t_with_factorization / CLOCKS_PER_SEC) << std::endl;
    std::cout << "\n" << std::endl;
}


TEST(PowerSeries3d, MixedHalfOpenAndClosedSegments) {
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
    double order_y = 18;
    double order_z = 18;
    double x = 10;
    double y = 15;
    double z = 14;
    auto myChebyshevSeries3dYZOpenedXClosed = ChebyshevSeries3dYZOpenedXClosed<double>(&myFunction3d, xmin, xmax, ymin,
                                                                                       zmin, order_x, order_y, order_z);
    myChebyshevSeries3dYZOpenedXClosed.Compute_aijk();

    // From Chebyshev series to power series.
    auto bijk = myChebyshevSeries3dYZOpenedXClosed.Compute_bijk();
    auto myPowerSeries3dYZOpenedXClosed = PowerSeries3dYZOpenedXClosed<double>(bijk, xmin, xmax, ymin, zmin);

    // Tests.
    PrintHeader("X closed, YZ open - Function.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries3dYZOpenedXClosed.Evaluate(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("X closed, YZ open - x derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dYZOpenedXClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("X closed, YZ open - y derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dYZOpenedXClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("X closed, YZ open - z derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_z(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dYZOpenedXClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dYZOpenedXClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    std::cout << "\n" << std::endl;
}


TEST(PowerSeries3d, MixedHalfOpenAndClosedSegmentsBis) {
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
    double order_z = 18;
    double x = 10;
    double y = 9;
    double z = 15;
    auto myChebyshevSeries3dZOpenedXYClosed = ChebyshevSeries3dZOpenedXYClosed<double>(&myFunction3d, xmin, xmax, ymin,
                                                                                       ymax, zmin, order_x, order_y,
                                                                                       order_z);
    myChebyshevSeries3dZOpenedXYClosed.Compute_aijk();

    // From Chebyshev series to power series.
    auto bijk = myChebyshevSeries3dZOpenedXYClosed.Compute_bijk();
    auto myPowerSeries3dZOpenedXYClosed = PowerSeries3dZOpenedXYClosed<double>(bijk, xmin, xmax, ymin, ymax, zmin);

    // Tests.
    PrintHeader("XY closed, Z open - Function.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    double ana = myFunction3d.Evaluate(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    double chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    double power = myPowerSeries3dZOpenedXYClosed.Evaluate(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("XY closed, Z open - x derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivative_x(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("XY closed, Z open - y derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivative_y(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("XY closed, Z open - z derivative.");
    std::cout << "Point: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_z(x, y, z);
    std::cout << "Analytical: " << ana << std::endl;
    chebyshev = myChebyshevSeries3dZOpenedXYClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Tripe Chebyshev series approximation: " << chebyshev << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((chebyshev - ana) / ana) << std::endl;
    power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivative_z(x, y, z);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, chebyshev, 10e-3));
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("XY closed, Z open - Function for z tending to +infinity.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction3d.Evaluate_zinf(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    power = myPowerSeries3dZOpenedXYClosed.Evaluate_zinf(x, y);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Relative error (%): " << 100 * abs((power - ana) / ana) << std::endl;
    EXPECT_TRUE(IsClose(ana, power, 10e-3));

    PrintHeader("XY closed, Z open - x derivative. for z tending to +infinity.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_x_zinf(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivative_x_zinf(x, y);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Absolute error (%): " << 100 * abs(power - ana) << std::endl; // Relative error infinite because ana = 0.
    EXPECT_TRUE(IsClose(abs(power - ana), 0., 10e-3, 10e-3));

    PrintHeader("XY closed, Z open - y derivative. for z tending to +infinity.");
    std::cout << "Point: (x, y) = (" << x << ", " << y << ")" << std::endl;
    ana = myFunction3d.Evaluate_derivative_y_zinf(x, y);
    std::cout << "Analytical: " << ana << std::endl;
    power = myPowerSeries3dZOpenedXYClosed.Evaluate_derivative_y_zinf(x, y);
    std::cout << "Triple Power series approximation: " << power << std::endl;
    std::cout << "Absolute error (%): " << 100 * abs(power - ana) << std::endl; // Relative error infinite because ana = 0.
    EXPECT_TRUE(IsClose(abs(power - ana), 0., 10e-3, 10e-3));
}
