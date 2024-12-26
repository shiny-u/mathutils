//
// Created by frongere on 25/01/18.
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


TEST(VectorN, OldMain) {
    // ========================================================================
    //    Initialization.
    // ========================================================================

    PrintHeader("Initialization");

    PrintInfo("To 0:");
    auto vec = VectorN<double>(10);
    std::cout << vec.transpose() << "\n\n";

    PrintInfo("To 1:");
    vec.SetOnes();
    std::cout << vec.transpose() << "\n\n";

    // ========================================================================
    //    Randomize.
    // ========================================================================

    PrintHeader("Randomize");
    vec.Randomize();
    std::cout << vec.transpose() << "\n\n";

    // ========================================================================
    //    Sorting.
    // ========================================================================

    PrintHeader("Sorting");

    PrintInfo("Sorting in ascending order:");
    vec.Sort();
    std::cout << vec.transpose() << "\n\n";

    PrintInfo("Sorting in descending order:");
    vec.Sort(false);
    std::cout << vec.transpose() << "\n\n";


    // ========================================================================
    //    Multiplying matrix by vector.
    // ========================================================================

    PrintHeader("Multiplying");

    PrintInfo("Multiplying by identity matrix:");
    auto matrix = MatrixMN<double>(10, 10);
    matrix.SetIdentity();

    VectorN<double> vec1 = matrix * vec;
    std::cout << vec1.transpose() << "\n\n";
    EXPECT_TRUE(vec.IsEqual(vec1));

    matrix.Randomize();
    VectorN<double> vec2 = matrix * vec;
//    vec1 = matrix * vec; // Does not work !!
    PrintInfo("Multiplying by a random matrix:");
    std::cout << vec2.transpose() << "\n\n";

    // Complex expression
    PrintInfo("Complex expression:");
    auto matrix2 = MatrixMN<double>(10, 10);
    matrix2.Randomize();
    VectorN<double> vec3 = matrix2 * vec + matrix * vec1 - vec2;
    std::cout << vec3.transpose() << "\n\n";

    // ========================================================================
    //    Generating a matrix from a vector.
    // ========================================================================

    PrintHeader("Generating a matrix from a vector");

    PrintInfo("Square matrix by function:");
    auto squareMat_function = GetMatrixSquare(vec3);
    std::cout << squareMat_function << "\n\n";

    PrintInfo("Square matrix by method:");
    auto squareMat_method = vec3.GetMatrixSquare();
    std::cout << squareMat_method << "\n\n";
    EXPECT_TRUE(squareMat_function.IsEqual(squareMat_method));

    PrintInfo("Diagonal to use:");
    std::cout << vec3.transpose() << "\n\n";

    PrintInfo("Diagonal matrix by function:");
    auto diagMat_function = GetDiagonalMatrix(vec3);
    std::cout << diagMat_function << "\n\n";

    PrintInfo("Diagonal matrix by method:");
    auto diagMat_method = vec3.GetDiagonalMatrix();
    std::cout << diagMat_method << "\n\n";
    EXPECT_TRUE(diagMat_function.IsEqual(diagMat_method));
}
