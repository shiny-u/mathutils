//
// Created by frongere on 29/11/17.
//

#include <fstream>
#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


#define N 101


using namespace mathutils;


TEST(Interp2d, OldMain) {
    // TODO Add more assertions.

    unsigned long it;

    // -----------------------------------------------
    // Test with one cell
    // -----------------------------------------------
    /**
    auto x = std::make_shared<std::vector<double>>();
    x->push_back(0.);
    x->push_back(1.);
    auto y = std::make_shared<std::vector<double>>();
    y->push_back(0.);
    y->push_back(1.);

    auto data = std::make_shared<std::vector<double>>();
    data->push_back(0.);
    data->push_back(1.);
    data->push_back(1.);
    data->push_back(2.);

    FrInterp2dLinear<double> interpolator;

    interpolator.Initialize(x, y, data);

    std::cout << " Eval : " << interpolator.Eval(0.75,0.5) << std::endl;
    */
    // ------------------------------------------------
    // Test with four cells
    // ------------------------------------------------
    /**
    auto x = std::make_shared<std::vector<double>>();
    x->push_back(0.);
    x->push_back(1.);
    x->push_back(2.);
    auto y = std::make_shared<std::vector<double>>();
    y->push_back(0.);
    y->push_back(1.);
    y->push_back(2.);

    auto data = std::make_shared<std::vector<double>>();
    data->push_back(0.);
    data->push_back(1.);
    data->push_back(0.);
    data->push_back(1.);
    data->push_back(2.);
    data->push_back(1.);
    data->push_back(2.);
    data->push_back(3.);
    data->push_back(2.);

    FrInterp2dLinear<double> interpolator;

    interpolator.Initialize(x, y, data);

    std::cout << " Eval : " << interpolator.Eval(2.,1.) << std::endl;

    auto x_interp = linspace(0., 2., 50);
    auto y_interp = linspace(0., 2., 50);
    */
    // ------------------------------------------------
    // Test with large amount of data and sinusoidal function
    // ------------------------------------------------

    auto x = std::make_shared<std::vector<double>>(linspace(M_PI, 4 * M_PI, N - 1));
    auto y = std::make_shared<std::vector<double>>(linspace(M_PI, 4 * M_PI, N - 1));

    auto data = std::make_shared<std::vector<double>>();

    double valx, valy;
    for (unsigned long i = 0; i < x->size(); i++) {
        valx = sin(x->at(i));
        for (unsigned long j = 0; j < y->size(); j++) {
            valy = sin(y->at(j));
            data->push_back(valx * valy);
        }
    }

    Interp2dLinear<double> interpolator;

    interpolator.Initialize(x, y, data);

    double yExpected = sin(5.333) * sin(5.000);
    auto y0 = interpolator.Eval(5.333, 5.000);

    std::cout << "" << std::endl;
    std::cout << "Test with large amount of data" << std::endl;
    std::cout << "Analytical : " << yExpected << std::endl;
    std::cout << "Eval : " << y0 << std::endl;

    EXPECT_NEAR(yExpected, y0, 1e-3);

    auto x_interp = linspace(M_PI, 4 * M_PI, 3 * (N - 1));
    auto y_interp = linspace(M_PI, 4 * M_PI, 3 * (N - 1));


    // ------------------------------------------------
    // Compute grid interpolation
    // ------------------------------------------------

    auto nx_interp = x_interp.size();
    auto ny_interp = y_interp.size();
    auto nx_ny = nx_interp * ny_interp;

    std::vector<std::vector<double>> vcoord;
    vcoord.reserve(x_interp.size() * y_interp.size());


    for (auto xi : x_interp) {
        for (auto yj : y_interp) {
            vcoord.push_back({xi, yj});
        }
    }

    auto tmp = vcoord.size();

    auto data_interp = interpolator.Eval(vcoord);

    // ----------------------------------------------
    // Save data in csv
    // ----------------------------------------------
    /**
    std::ofstream myfile;
    myfile.open("data.csv");

    for (unsigned long i = 0; i < x->size(); i++) {
        for (unsigned long j = 0; j < y->size(); j++) {
            it = i * y->size() + j;
            myfile << x->at(i) << ";" << y->at(j) << ";" << data->at(it) << std::endl;
        }
    }

    myfile.close();
    */
    // -----------------------------------------------
    // Save data_interp in csv
    // -----------------------------------------------

    std::ofstream myfile2;
    myfile2.open("data_interp.csv");

    for (unsigned long i = 0; i < x_interp.size(); i++) {
        for (unsigned long j = 0; j < y_interp.size(); j++) {
            it = i * y_interp.size() + j;
            myfile2 << x_interp.at(i) << ";" << y_interp.at(j) << ";" << data_interp.at(it) << std::endl;
        }
    }

    myfile2.close();
}
