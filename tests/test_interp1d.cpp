//
// Created by frongere on 16/11/17.
//

#include <iostream>
#include "MathUtils/MathUtils.h"
#include "gtest/gtest.h"

#include <gtest/gtest.h>


#define N 100


using namespace mathutils;


TEST(Interp1d, Double) {

  // Building the x coords as a shared pointer
  auto x = std::make_shared<std::vector<double>>(
      linspace(M_PI_2, 4 * M_PI, N + 1)
  );

  // Building the data
  auto y = std::make_shared<std::vector<double>>();
  y->reserve(x->size());
  double val;
  for (unsigned long i = 0; i < x->size(); i++) {
    val = sin(x->at(i));
    y->push_back(val);
  }

  // Create the linear interpolator
  Interp1dLinear<double, double> interpolator;
  interpolator.Initialize(x, y);

  // Test of the Eval method for one scalar
  EXPECT_EQ(interpolator.Eval(5.3333), interpolator(5.3333));
  EXPECT_NEAR(sin(5.3333), interpolator(5.3333), 1E-2);

  // Test exit if outside of ranges
  EXPECT_EXIT(interpolator(0.), ::testing::ExitedWithCode(1), ".*");

  // Test for a vector of x coords
  auto x_interp = linspace(M_PI_2, 4 * M_PI, 1000 * N);
  // Using only the overloaded call operator for vector values
  auto y_interp = interpolator(x_interp);

  // Interpolator saturates outside the ranges
  Interp1dLinearSaturate<double, double> saturator;
  saturator.Initialize(x, y);

  EXPECT_EQ(saturator(MU_PI_2), saturator(0.));
  EXPECT_EQ(saturator(4 * MU_PI), saturator(5 * MU_PI));

  // Linear extrapolator
  Interp1dLinearExtrapolate<double, double> extrapolator;
  extrapolator.Initialize(x, y);

  // test below x0
  auto a = (sin(x->at(1)) - sin(x->at(0))) / (x->at(1) - x->at(0));
  auto b = sin(x->at(0)) - a * x->at(0);
  EXPECT_NEAR(b, extrapolator(0.), 1E-8);

  // test above xn
  a = (sin(x->at(N)) - sin(x->at(N - 1))) / (x->at(N) - x->at(N - 1));
  b = sin(x->at(N)) - a * (x->at(N) - 5 * MU_PI);
  EXPECT_NEAR(b, extrapolator(5 * MU_PI), 1E-8);

}


TEST(Interp1d, Complex) {

  // Building the x coords as a shared pointer
  auto x = std::make_shared<std::vector<double>>(
      linspace(M_PI, 4 * M_PI, N + 1)
  );

  // Building the data
  auto y = std::make_shared<std::vector<std::complex<double>>>();
  y->reserve(x->size());
  std::complex<double> val;
  for (unsigned long i = 0; i < x->size(); i++) {
    val = exp(MU_JJ * x->at(i));
    y->push_back(val);
  }

  // Create the interpolation
  Interp1dLinear<double, std::complex<double>> interpolator;

  interpolator.Initialize(x, y);

  // Test of the Eval method for one scalar
  EXPECT_EQ(interpolator.Eval(5.3333), interpolator(5.3333));
  EXPECT_NEAR(exp(MU_JJ * 5.3333).real(), interpolator(5.3333).real(), 1E-2);
  EXPECT_NEAR(exp(MU_JJ * 5.3333).imag(), interpolator(5.3333).imag(), 1E-2);

  // Test for a vector of x coords
  auto x_interp = linspace(M_PI, 4 * M_PI, 1000 * N);
  // Using only the overloaded call operator for vector values
  auto y_interp = interpolator(x_interp);
}


TEST(Interp1d, Grad) {

    // Building the x coords as a shared pointer
    std::vector<double> x{0.0, 10.0, 100.0, 1000.0};
    auto x_coord = std::make_shared<std::vector<double>>(x);

    // Building the data
    std::vector<double> y{10.0, 10.0, 100.0, 2000.0};
    auto y_coord = std::make_shared<std::vector<double>>(y);

    // Create the interpolation
    Interp1dLinear<double, double> interpolator;

    interpolator.Initialize(x_coord, y_coord);


    EXPECT_NEAR(interpolator.GradEval(0.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(1.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(9.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(50.0),1.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(500.0),1900.0/900.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(1000.0),1900.0/900.0,1e-5);
}

TEST(Interp1dSaturate, Grad) {

    // Building the x coords as a shared pointer
    std::vector<double> x{0.0, 10.0, 100.0, 1000.0};
    auto x_coord = std::make_shared<std::vector<double>>(x);

    // Building the data
    std::vector<double> y{10.0, 10.0, 100.0, 2000.0};
    auto y_coord = std::make_shared<std::vector<double>>(y);

    // Create the interpolation
    Interp1dLinearSaturate<double, double> interpolator;

    interpolator.Initialize(x_coord, y_coord);


    EXPECT_NEAR(interpolator.GradEval(-100.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(0.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(1.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(9.0),0.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(50.0),1.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(500.0),1900.0/900.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(1000.0),1900.0/900.0,1e-5);
    EXPECT_NEAR(interpolator.GradEval(10000.0),0.0,1e-5);
}


TEST(Interp1dSaturate, ComplexGrad) {

    // Building the x coords as a shared pointer
    auto x = std::make_shared<std::vector<double>>(
            linspace(M_PI, 4 * M_PI, N + 1)
    );

    // Building the data
    auto y = std::make_shared<std::vector<std::complex<double>>>();
    y->reserve(x->size());
    std::complex<double> val;
    for (unsigned long i = 0; i < x->size(); i++) {
        val = exp(MU_JJ * x->at(i));
        y->push_back(val);
    }

    // Create the interpolation
    Interp1dLinearSaturate<double, std::complex<double>> interpolator;

    interpolator.Initialize(x, y);
    auto output = interpolator.GradEval(0);

    ASSERT_NEAR(output.real(),0.0,1e-5);
}

