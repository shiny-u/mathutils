//
// Created by bpenin on 06/06/19.
//

#include <iostream>
#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace mathutils;


TEST(Geometry, QuaternionWithEuler) {
    Eigen::Vector3d euler_angles(30*MU_PI/180,45*MU_PI/180,-60*MU_PI/180);

    EXPECT_NEAR(quaternionToEulerAnglesZYX(eulerAnglesZYXToQuaternion(euler_angles))(0), euler_angles(0), 1e-3);
    EXPECT_NEAR(quaternionToEulerAnglesZYX(eulerAnglesZYXToQuaternion(euler_angles))(1), euler_angles(1), 1e-3);
    EXPECT_NEAR(quaternionToEulerAnglesZYX(eulerAnglesZYXToQuaternion(euler_angles))(2), euler_angles(2), 1e-3);
}


TEST(Geometry, RotationMatrixWithEuler) {
    Eigen::Vector3d euler_angles(30*MU_PI/180,45*MU_PI/180,-60*MU_PI/180);

    EXPECT_NEAR(rotationMatrixToEulerAnglesZYX(eulerAnglesZYXToRotationMatrix(euler_angles))(0), euler_angles(0), 1e-3);
    EXPECT_NEAR(rotationMatrixToEulerAnglesZYX(eulerAnglesZYXToRotationMatrix(euler_angles))(1), euler_angles(1), 1e-3);
    EXPECT_NEAR(rotationMatrixToEulerAnglesZYX(eulerAnglesZYXToRotationMatrix(euler_angles))(2), euler_angles(2), 1e-3);
}


TEST(Geometry, RotationMatrixWithQuaternion) {
    Eigen::Quaterniond quat(0.7,0,0,0.7);

    EXPECT_NEAR(RotationMatrixToQuaternion(quaternionToRotationMatrix(quat)).x(), quat.x(), 1e-3);
    EXPECT_NEAR(RotationMatrixToQuaternion(quaternionToRotationMatrix(quat)).y(), quat.y(), 1e-3);
    EXPECT_NEAR(RotationMatrixToQuaternion(quaternionToRotationMatrix(quat)).z(), quat.z(), 4e-3);
    EXPECT_NEAR(RotationMatrixToQuaternion(quaternionToRotationMatrix(quat)).w(), quat.w(), 4e-3);
}
