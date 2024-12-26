//
// Created by frongere on 16/11/17.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace mathutils;


TEST(Unit, OldMain) {

    EXPECT_NEAR(180. * DEG2RAD, M_PI, 1e-8);
    EXPECT_NEAR(M_PI * RAD2DEG, 180, 1e-8);

    // Pas tres utilse de tout faire...
}
