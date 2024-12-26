//
// Created by frongere on 16/11/17.
//

#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace mathutils;


TEST(Angles, Normalize)
{
    assert(Normalize__PI_PI(3*MU_PI_2) == -MU_PI_2);
    assert(Normalize__PI_PI(7*MU_PI_2) == -MU_PI_2);
    assert(Normalize__PI_PI(-3*MU_PI_2) == MU_PI_2);
    assert(Normalize__PI_PI(-7*MU_PI_2) == MU_PI_2);
    assert(Normalize__PI_PI(MU_2PI) == 0.);
    assert(Normalize__PI_PI(-M_PI) == M_PI);
    assert(Normalize__PI_PI(M_PI) == M_PI);

    assert(Normalize__180_180(3*90) == -90);
    assert(Normalize__180_180(7*90) == -90);
    assert(Normalize__180_180(-3*90) == 90);
    assert(Normalize__180_180(-7*90) == 90);
    assert(Normalize__180_180(360) == 0);
    assert(Normalize__180_180(-180) == 180);
    assert(Normalize__180_180(180) == 180);

    assert(Normalize_0_2PI(3*M_PI) == M_PI);
    assert(Normalize_0_2PI(-M_PI) == M_PI);
    assert(Normalize_0_2PI(0) == 0);
    assert(Normalize_0_2PI(MU_2PI) == 0);

    assert(Normalize_0_360(540) == 180);
    assert(Normalize_0_360(-180) == 180);
    assert(Normalize_0_360(0) == 0);
    assert(Normalize_0_360(360) == 0);
}
