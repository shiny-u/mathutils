//
// Created by pierre-yves on 14/01/20.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace std;
using namespace mathutils;


TEST(Integration2d, OldMain) {
    // This test checks the integration of a function over a triangle.
    // This function returns a VectorN because it is more difficult when the size is not fixed.

    // Vertices of a triangle.
    Vector3d<double> vertex_1 = Vector3d<double>(0., 0., 0.);
    Vector3d<double> vertex_2 = Vector3d<double>(1., 0.5, 0.);
    Vector3d<double> vertex_3 = Vector3d<double>(0., 1., 0.);

    // Order of the integration.
    int order = 2;

    // Definition of a function to integration.
    class IntegrandTest : public Integrand<VectorN<double>> {
    public:
        VectorN<double> Evaluate(const Vector3d<double> &x) const override {
            VectorN<double> value = VectorN<double>::Zero(1);
            value[0] = 2 - x(0) - 2 * x(1);
            return value; // f(x,y,z) = 2 - x - 2y.
        }
    };

    // Numerical integration.
    IntegrandTest myFunction;
    auto myIntegrator = Integration2dTriangle<VectorN<double>>(&myFunction, order);
    VectorN<double> numerical_result = myIntegrator.Compute(vertex_1, vertex_2, vertex_3);

    // Analytical integration.
    double analytical_result = 1./3.;

    std::cout << "" << std::endl;
    std::cout << "Analytical result: " << analytical_result << std::endl;
    std::cout << "Numerical result: " << numerical_result[0] << std::endl;
    std::cout << "Relative error: " << (100*(numerical_result[0] - analytical_result) / analytical_result) << std::endl;

    EXPECT_NEAR(analytical_result, numerical_result[0], 1e-13);
}
