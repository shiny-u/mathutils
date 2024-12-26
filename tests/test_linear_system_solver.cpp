//
// Created by pierre-yves on 22/10/19.
//

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


TEST(LinearSystemSolver, OldMain) {
    PrintHeader("Matrix66");

    // Matrix66.
    mathutils::Matrix66<double> A66;
    A66 << 2,5,7,9,1,-6,
            3,7,-4,8,-6,1,
            -10,13,24,-14,3,0,
            3,0,4,8,9,-6,
            1,7,9,3,-20,4,
            0,0,3,4,8,9;

    PrintInfo("The matrix A");
    std::cout << A66 << "\n\n";

    mathutils::Vector6d<double> b6;
    b6 << 1, 2, 3, 4, 5, 6;

    PrintInfo("The right-hand side vector b");
    std::cout << b6 << "\n\n";

    auto x6 = A66.LUSolver(b6);

    PrintInfo("The solution x from a LU decomposition");
    std::cout << x6 << "\n\n";

    PrintInfo("Verification (A * x - b = 0)");
    std::cout << A66 * x6 - b6 << "\n\n";
    Vector6d<double> Ax66 = A66 * x6;
    EXPECT_TRUE(Ax66.IsEqual(b6));

    // Matrix33.

    PrintHeader("Matrix33");

    mathutils::Matrix33<double> A33;
    A33 << 2,5,7,
            3,7,-4,
            -10,13,24;

    PrintInfo("The matrix A");
    std::cout << A33 << "\n\n";

    mathutils::Vector3d<double> b3;
    b3 << 1, 2, 3;

    PrintInfo("The right-hand side vector b");
    std::cout << b3 << "\n\n";

    auto x3 = A33.LUSolver(b3);

    PrintInfo("The solution x from a LU decomposition");
    std::cout << x3 << "\n\n";

    PrintInfo("Verification (A * x - b = 0)");
    std::cout << A33 * x3 - b3 << "\n\n";
    Vector3d<double> Ax33 = A33 * x3;
    EXPECT_TRUE(Ax33.IsEqual(b3));

    // MatrixMN.

    PrintHeader("MatrixMN");

    mathutils::MatrixMN<double> AMN = mathutils::MatrixMN<double>::Zero(6, 6);
    AMN = A66;

    PrintInfo("The matrix A");
    std::cout << A66 << "\n\n";

    mathutils::MatrixMN<double> bN = mathutils::MatrixMN<double>::Zero(6, 2);
    bN << 1, 1, 2, 2, 3, 3,
            4, 4, 5, 5, 6, 6;

    PrintInfo("The right-hand side vector b");
    std::cout << bN << "\n\n";

    auto xN = AMN.LUSolver(bN);

    PrintInfo("The solution x from a LU decomposition");
    std::cout << xN << "\n\n";

    PrintInfo("Verification (A * x - b = 0)");
    std::cout << AMN * xN - bN << "\n\n";
    mathutils::MatrixMN<double> AxMN = AMN * xN;
    EXPECT_TRUE(AxMN.IsEqual(bN));
}
