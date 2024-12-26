//
// Created by pierre-yves on 31/10/20.
//

#include <iostream>
#include "MathUtils/MathUtilsBoost.h"
#include "gtest/gtest.h"

using namespace mathutils;

void PrintHeader(std::string title) {
  std::cout << "\n=====================================================================" << std::endl;
  std::cout << "    " << title << std::endl;
  std::cout << "=====================================================================" << std::endl;
}

void PrintInfo(std::string info) {
  std::cout << info << ":" << std::endl;
}

// This test checks the computation of the special functions.

// ========================================================================
//                          Bessel functions
// ========================================================================

// Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_j
TEST(BoostFunctions, Cyl_Bessel_first_kind) {
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0 and x = 1.2345");
  double value_J0 = Cyl_Bessel_first_kind<int, double, double>(0, 1.2345);
  std::cout << "Method Cyl_Bessel_first_kind: " << value_J0 << std::endl;
  std::cout << "Exact value: 0.653792" << std::endl;
  EXPECT_NEAR(value_J0, 0.653792, 1e-6);
}

TEST(BoostFunctions, Cyl_Bessel_second_kind) {
  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_neumann
//  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0.5 and x = 1.2345");
  double value_Y05 = Cyl_Bessel_second_kind<double, double, double>(0.5, 1.2345);
  std::cout << "Method Cyl_Bessel_second_kind: " << value_Y05 << std::endl;
  std::cout << "Exact value: -0.236973" << std::endl;
  EXPECT_NEAR(value_Y05, -0.236973, 1e-6);
}

TEST(BoostFunctions, Cyl_modified_Bessel_first_kind) {
  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_i
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0 and x = 1.2345");
  double value_I0 = Cyl_modified_Bessel_first_kind<double, double, double>(0, 1.2345);
  std::cout << "Method Cyl_modified_Bessel_first_kind: " << value_I0 << std::endl;
  std::cout << "Exact value: 1.41886" << std::endl;
  EXPECT_NEAR(value_I0, 1.41886, 1e-6);
}

TEST(BoostFunctions, Cyl_modified_Bessel_second_kind) {
  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_k
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0.5 and x = 1.2345");
  double value_K05 = Cyl_modified_Bessel_second_kind<double, double, double>(0.5, 1.2345);
  std::cout << "Method Cyl_modified_Bessel_second_kind: " << value_K05 << std::endl;
  std::cout << "Exact value: 0.32823" << std::endl;
  EXPECT_NEAR(value_K05, 0.32823, 1e-6);
}

TEST(BoostFunctions, Cyl_Hankel_first_kind) {
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 0 and x = 1.2345");
  std::complex<double> value_H0 = Cyl_Hankel_first_kind<int, double>(0, 1.2345);
  std::cout << "Method Cyl_Hankel_first_kind: " << value_H0 << std::endl;
  std::cout << "Exact value: 0.653792 + 0.249072i" << std::endl;
  EXPECT_NEAR(value_H0.real(), 0.653792, 1e-6);
  EXPECT_NEAR(value_H0.imag(), 0.249072, 1e-6);
}

// ========================================================================
//                       Legendre polynomials
// ========================================================================

TEST(BoostFunctions, Legendre_polynomial) {
  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/legendre
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 3 and x = 0.25");
  double value_Legendre_P3 = Legendre_polynomial<int, double>(3, 0.25);
  std::cout << "Method Legendre_polynomial: " << value_Legendre_P3 << std::endl;
  std::cout << "Exact value: " << 0.5 * (5 * 0.25 * 0.25 * 0.25 - 3 * 0.25) << std::endl; // P_3(x) = 1/2*(5x^3 - 3x)
  EXPECT_NEAR(value_Legendre_P3, 0.5 * (5 * 0.25 * 0.25 * 0.25 - 3 * 0.25), 1e-10);
}
// ========================================================================
//                 Derivative of the Legendre polynomials
// ========================================================================

TEST(BoostFunctions, Legendre_polynomial_derivative) {
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 3 and x = 0.25");
  double value_derivative_Legendre_P3 = Legendre_polynomial_derivative<int, double>(3, 0.25);
  std::cout << "Method Legendre_derivative: " << value_derivative_Legendre_P3 << std::endl;
  std::cout << "Exact value: " << 0.5 * (15 * 0.25 * 0.25 - 3) << std::endl; // 1/2*(15x^2 - 3)
  EXPECT_NEAR(value_derivative_Legendre_P3, 0.5 * (15 * 0.25 * 0.25 - 3), 1e-10);
}
// ========================================================================
//                            Factorial
// ========================================================================

TEST(BoostFunctions, Factorial) {
//  PrintInfo("Comparison of the value returned by the method and an exact value for n = 9");
  double value_factorial = Factorial<int, double>(9);
  std::cout << "Method Factorial: " << value_factorial << std::endl;
  std::cout << "Exact value: " << 362880. << std::endl; // 9! = 362880
  EXPECT_NEAR(value_factorial, 362880., 1e-10);
}
// ========================================================================
//                        Exponential integral Ei
// ========================================================================

TEST(BoostFunctions, Exponential_integral) {
  // Source: https://en.cppreference.com/w/cpp/numeric/special_functions/expint
  PrintInfo("Comparison of the value returned by the direct method and an exact value for x = 1");
  double value_Ei = Ei<double>(1.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Exact value: 1.89511781" << std::endl;
  EXPECT_NEAR(value_Ei, 1.89511781, 1e-5);

  double value_Ei_approx = Ei_approximation<double>(1.);
  std::cout << "Method Ei_approximation: " << value_Ei_approx << std::endl;
  std::cout << "Exact value: 1.89511781" << std::endl;
  EXPECT_NEAR(value_Ei_approx, 1.89511781, 1e-5);

  PrintInfo("\nComparison of the value returned by the two methods for evaluating Ei");
  PrintInfo("x = 3");
  value_Ei = Ei<double>(3.);
  value_Ei_approx = Ei_approximation<double>(3.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  EXPECT_NEAR(abs(value_Ei - value_Ei_approx), 0., 1e-10);

  PrintInfo("\nx = 9");
  value_Ei = Ei<double>(9.);
  value_Ei_approx = Ei_approximation<double>(9.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  EXPECT_NEAR(abs(value_Ei - value_Ei_approx), 0., 1e-9);

  PrintInfo("\nx = 18");
  value_Ei = Ei<double>(18.);
  value_Ei_approx = Ei_approximation<double>(18.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  EXPECT_NEAR(abs(value_Ei - value_Ei_approx), 0., 1e-10);

  PrintInfo("\nx = 30");
  value_Ei = Ei<double>(30.);
  value_Ei_approx = Ei_approximation<double>(30.);
  std::cout << "Method Ei: " << value_Ei << std::endl;
  std::cout << "Method Ei_approximation: " << value_Ei_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_Ei - value_Ei_approx) << std::endl;
  EXPECT_NEAR(abs(value_Ei - value_Ei_approx), 0., 1e-10);
}

TEST(BoostFunctions, Exponential_integral_direct) {
  PrintInfo(
      "\nComparison of the value returned by the method for evaluating exp(-x)Ei(x) compared to a direct computation.");
  PrintInfo("x = 3");
  double value_expEi = exp(-3.) * Ei<double>(3.);
  double value_expEi_approx = expEi<double>(3.);
  std::cout << "Direct: " << value_expEi << std::endl;
  std::cout << "Method expEi: " << value_expEi_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_expEi - value_expEi_approx) << std::endl;
  EXPECT_NEAR(abs(value_expEi - value_expEi_approx), 0., 1e-10);

  PrintInfo("x = 25");
  value_expEi = exp(-25.) * Ei<double>(25.);
  value_expEi_approx = expEi<double>(25.);
  std::cout << "Direct: " << value_expEi << std::endl;
  std::cout << "Method expEi: " << value_expEi_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_expEi - value_expEi_approx) << std::endl;
  EXPECT_NEAR(abs(value_expEi - value_expEi_approx), 0., 1e-10);

  PrintInfo("x = 50");
  value_expEi = exp(-50.) * Ei<double>(50.);
  value_expEi_approx = expEi<double>(50.);
  std::cout << "Direct: " << value_expEi << std::endl;
  std::cout << "Method expEi: " << value_expEi_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_expEi - value_expEi_approx) << std::endl;
  EXPECT_NEAR(abs(value_expEi - value_expEi_approx), 0., 1e-10);

  PrintInfo("x = 2000");
  value_expEi_approx = expEi<double>(2000.);
  std::cout << "Direct: Numerical overflow" << std::endl;
  std::cout << "Method expEi: " << value_expEi_approx << std::endl;
}

// ========================================================================
//                        Exponential integral E1
// ========================================================================

TEST(BoostFunctions, Exponential_integral_E1) {
  // Source: 1964_ABRAMOWITZ_Handbook_of_mathematical_functions_with_formulas_graphs_and_mathematical_tables (p239).
  PrintInfo("Comparison of the value returned by the direct method and an exact value for x = 1");
  double value_E1 = E1<double>(1.);
  std::cout << "Method E1: " << value_E1 << std::endl;
  std::cout << "Exact value: 0.219383934" << std::endl;
  EXPECT_NEAR(value_E1, 0.219383934, 1e-5);

  double value_E1_approx = E1_approximation<double>(1.);
  std::cout << "Method E1_approximation: " << value_E1_approx << std::endl;
  std::cout << "Exact value: 0.219383934" << std::endl;
  EXPECT_NEAR(value_E1_approx, 0.219383934, 1e-5);

  PrintInfo("\nComparison of the value returned by the two methods for evaluating Ei");
  PrintInfo("x = 0.5");
  value_E1 = E1<double>(0.5);
  value_E1_approx = E1_approximation<double>(0.5);
  std::cout << "Method E1: " << value_E1 << std::endl;
  std::cout << "Method E1_approximation: " << value_E1_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_E1 - value_E1_approx) << std::endl;
  EXPECT_NEAR(abs(value_E1 - value_E1_approx), 0., 1e-8);

  PrintInfo("x = 3");
  value_E1 = E1<double>(3.);
  value_E1_approx = E1_approximation<double>(3.);
  std::cout << "Method E1: " << value_E1 << std::endl;
  std::cout << "Method E1_approximation: " << value_E1_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_E1 - value_E1_approx) << std::endl;
  EXPECT_NEAR(abs(value_E1 - value_E1_approx), 0., 1e-8);

  PrintInfo("x = 20");
  value_E1 = E1<double>(20.);
  value_E1_approx = E1_approximation<double>(20.);
  std::cout << "Method E1: " << value_E1 << std::endl;
  std::cout << "Method E1_approximation: " << value_E1_approx << std::endl;
  std::cout << "|Delta| = " << abs(value_E1 - value_E1_approx) << std::endl;
  EXPECT_NEAR(abs(value_E1 - value_E1_approx), 0., 1e-8);

}

// ========================================================================
//                       Chebyshev polynomials
// ========================================================================

TEST(BoostFunctions, Chebyshev_polynomial) {
  // Source: https://fr.wikipedia.org/wiki/Polyn%C3%B4me_de_Tchebychev
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 5 and x = 0.25");
  double value_Chebyshev_T5_1 = Chebyshev_polynomial<int, double>(5, 0.25);
  std::cout << "Method Chebyshev_polynomial: " << value_Chebyshev_T5_1 << std::endl;
  double value_Chebyshev_T3 = Chebyshev_polynomial<int, double>(3, 0.25);
  double value_Chebyshev_T4 = Chebyshev_polynomial<int, double>(4, 0.25);
  double value_Chebyshev_T5_2 = Chebyshev_polynomial_next<double>(0.25, value_Chebyshev_T4, value_Chebyshev_T3);
  std::cout << "Method Chebyshev_polynomial_next: " << value_Chebyshev_T5_2 << std::endl;
  std::cout << "Exact value: " << 16. * pow(0.25, 5.) - 20. * pow(0.25, 3.) + 5. * 0.25
            << std::endl; // T_5(x) = 16x^5 - 20x^3 + 5x.
  EXPECT_NEAR(value_Chebyshev_T5_1, 16. * pow(0.25, 5.) - 20. * pow(0.25, 3.) + 5. * 0.25, 1e-10);
  EXPECT_NEAR(value_Chebyshev_T5_2, 16. * pow(0.25, 5.) - 20. * pow(0.25, 3.) + 5. * 0.25, 1e-10);
}
// ========================================================================
//                   Derivate of Chebyshev polynomials
// ========================================================================

TEST(BoostFunctions, Chebyshev_polynomial_derivative) {
  PrintInfo("Comparison of the value returned by the method and an exact value for n = 5 and x = 0.25");
  double value_Chebyshev_derivative_T5 = Chebyshev_polynomial_derivative<int, double>(5, 0.25);
  std::cout << "Method Chebyshev_polynomial_derivative: " << value_Chebyshev_derivative_T5 << std::endl;
  std::cout << "Exact value: " << 80. * pow(0.25, 4.) - 60. * pow(0.25, 2.) + 5
            << std::endl; // T'_5(x) = 80x^4 - 60x^2 + 5.
  EXPECT_NEAR(value_Chebyshev_derivative_T5, 80. * pow(0.25, 4.) - 60. * pow(0.25, 2.) + 5, 1e-10);
}
