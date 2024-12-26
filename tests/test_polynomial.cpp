//
// Created by pierre-yves on 03/03/2021.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

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


TEST(Polynomial, OldMain) {
    // This test checks the application of the Horner's method for the computation of a polynomial.

    std::vector<double> coefficient = {2.3, 3.4, -5, 7.8, 10, 30};
    double x0 = 3.;

    PrintHeader("Horder's method on a polynomial");
    double polynomial = 0.;
    double xn = 1.;
    for (int n = 0; n < coefficient.size(); ++n){
        polynomial += xn * coefficient.at(n);
        xn *= x0;
    }
    std::cout << "Iterative method: " << polynomial << std::endl;
    double horner = Horner<double>(coefficient, x0);
    std::cout << "Horner's method: " << horner << std::endl;
    EXPECT_NEAR(polynomial, horner, 1e-8);

    PrintHeader("Horder's method on the derivative of a polynomial");
    double derivate = 0.;
    xn = 1.;
    double xn_minus = 0.;
    for (int n = 0; n < coefficient.size(); ++n){
        derivate += n * xn_minus * coefficient.at(n);
        xn_minus = xn;
        xn *= x0;
    }
    std::cout << "Iterative method: " << derivate << std::endl;
    double horner_derivative = Horner_derivative<double>(coefficient, x0);
    std::cout << "Horner's method: " << horner_derivative << std::endl;
    EXPECT_NEAR(derivate, horner_derivative, 1e-8);
}
