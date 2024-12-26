//
// Created by frongere on 16/11/17.
//


//#include <iostream>
#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace mathutils;


TEST(LookupTable, OldMain) {
    auto x = linspace<double>(0, 100, 100);
    auto cx = linspace<double>(10, 90, 100);
    auto cy = linspace<double>(0.2, 87., 100);

    auto lookup = LookupTable1D<double>(LINEAR);
    lookup.SetX(x);
    lookup.AddY("cx", cx);
    lookup.AddY("cy", cy);

    auto lookup_2 = LookupTable1D_DataContained<double>(LINEAR);
    lookup_2.SetX(x);
    lookup_2.AddY("cx", cx);
    lookup_2.AddY("cy", cy);

    auto xInterp = linspace<double>(0.1, 99.9, 200);

    lookup.Eval("cx", xInterp);
    auto lookup_3(lookup_2);

    auto lookup_4 = LookupTable1D_DataContained<double>(LINEAR_EXTRAPOLATE);

    lookup_4 = lookup_2;

    ASSERT_TRUE(lookup_4.get_interp_method() == LINEAR);

    for (auto xi : xInterp) {

        ASSERT_NEAR(10 + 80 * xi/100, lookup.Eval("cx", xi), 1e-8);
        ASSERT_NEAR(0.2 + 86.8 * xi/100, lookup.Eval("cy", xi), 1e-8);

        ASSERT_NEAR(10 + 80 * xi/100, lookup_2.Eval("cx", xi), 1e-8);
        ASSERT_NEAR(0.2 + 86.8 * xi/100, lookup_2.Eval("cy", xi), 1e-8);

        ASSERT_NEAR(10 + 80 * xi/100, lookup_3.Eval("cx", xi), 1e-8);
        ASSERT_NEAR(0.2 + 86.8 * xi/100, lookup_3.Eval("cy", xi), 1e-8);
    }

}
