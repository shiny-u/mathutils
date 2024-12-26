//
// Created by frongere on 14/12/17.
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


TEST(MatrixMN, OldMain) {
    // ========================================================================
    //    Matrix creation, initialization and extractions (indexing)
    // ========================================================================
    PrintHeader("Instantiating a matrix with 6 rows and 3 columns");

    // Instatiating a matrix with its dimensions (it is set to zeros by default)
    auto myMatrix = MatrixMN<double>(6, 3);

    // Displaying the matrix
    PrintInfo("Initial Matrix");
    std::cout << myMatrix << "\n\n";

    // Getting dimensions of the matrix
    std::cout << "Nb Rows : " << myMatrix.GetNbRows() << std::endl << std::endl;
    std::cout << "Nb Cols : " << myMatrix.GetNbCols() << std::endl << std::endl;

    PrintInfo("Randomizing coefficients");
    // Randomizing the coefficients
    myMatrix.Randomize();
    std::cout << myMatrix << "\n\n";

    // Getting individual coefficients (different means)
    PrintInfo("Getting and setting some individual coefficients by indices");
    std::cout << "Get Coeff 0, 0 : " << myMatrix.at(0, 0) << "\n\n";
    std::cout << "Get Coeff 0, 2 : " << myMatrix.at(0, 2) << "\n\n";
    std::cout << "Get Coeff 0, 0 : " << myMatrix(0, 0) << "\n\n";

    myMatrix.at(0, 0) = 0;
    std::cout << "Set Coeff 0, 0 (0): " << myMatrix(0, 0) << "\n\n";

    myMatrix(1, 1) = 11;
    std::cout << "Set Coeff 1, 1 (11): " << myMatrix(1, 1) << "\n\n";

    myMatrix(0, 1) = 1;
    std::cout << "Set Coeff 0, 1 (1): " << myMatrix(0, 1) << "\n\n";

    myMatrix(1, 2) = 12;
    std::cout << "Set Coeff 1, 2 (12): " << myMatrix(1, 2) << "\n\n";

    myMatrix(3, 0) = 30;
    std::cout << "Set Coeff 3, 0 (30): " << myMatrix(3, 0) << "\n\n";

    myMatrix.at(4, 0) = 40;
    std::cout << "Set Coeff 4, 0 (40): " << myMatrix(4, 0) << "\n\n";

    std::cout << myMatrix << "\n\n";

    PrintInfo("Matrix transposition");
    // Transpose the matrix inplace
    myMatrix.Transpose();
    std::cout << myMatrix << "\n\n";

    std::cout << "Nb Rows : " << myMatrix.GetNbRows() << std::endl << std::endl;
    std::cout << "Nb Cols : " << myMatrix.GetNbCols() << std::endl << std::endl;

    PrintInfo("Diagonal extraction");
    // Get the diagonal
    auto diag = myMatrix.GetDiag();
    std::cout << diag << "\n\n";

    PrintInfo("Column extraction by index");
    // Get a column
    auto col2 = myMatrix.GetColumn(2);
    std::cout << col2 << "\n\n";

    PrintInfo("Row extraction by index");
    // Get a row
    auto row1 = myMatrix.GetRow(1);
    std::cout << row1 << "\n\n";

    PrintInfo("Setting coeffs to 0");
    // Setting every coeffs to 0
    myMatrix.SetNull();
    std::cout << myMatrix << "\n\n";

    PrintInfo("Set identity matrix");
    // Making it identity (even in square matrices)
    myMatrix.SetIdentity();
    std::cout << myMatrix << "\n\n";

    // ========================================================================
    //    Matrix multiplication
    // ========================================================================
    PrintHeader("Multiplying matrix by matrix");
    auto matrix2 = MatrixMN<double>(10, 10);
    matrix2.Randomize();
    auto matrix = MatrixMN<double>(10, 10);
    matrix.Randomize();
    MatrixMN<double> matrix3 = matrix * matrix2;
    std::cout << matrix3 << "\n\n";

    // ========================================================================
    //    Square matrix inversion
    // ========================================================================

    PrintHeader("Inversion of square matrices");
    MatrixMN<double> matSquare(6, 6);
    matSquare.Randomize();
    PrintInfo("The square Matrix to inverse");
    std::cout << matSquare << "\n\n";

    PrintInfo("Inverse of this matrix inplace");
    // Inverse of the matrix (inplace)
    auto initialMatrix = matSquare;
    matSquare.Inverse();
    std::cout << matSquare << "\n\n";

    PrintInfo("Inverse of the same matrix by copy (we come back to the original...)");
    // Inverse of the Matrix (copy)
    auto invMatrix = matSquare.GetInverse();
    std::cout << invMatrix << "\n\n"; // Must get back the initial matrix

    PrintInfo("Verifications ( (A^-1)^-1 - A = 0)");
    std::cout << invMatrix - initialMatrix << "\n\n"; // Must be zero
    EXPECT_TRUE(invMatrix.IsEqual(initialMatrix));

    PrintInfo("Verifications (A^-1 * A = I)");
    MatrixMN<double> multInvInit = matSquare * initialMatrix;
    std::cout << multInvInit << "\n\n";   // Must be identity
    EXPECT_TRUE(multInvInit.IsIdentity(1e-8));

//    // Is the matrix positive semi definite
//    auto isPSD = matSquare.IsPositiveSemiDefinite();


    // ========================================================================
    //    QR Decomposition
    // ========================================================================

    PrintHeader("Decompositions");
    PrintHeader("Thin QR Decomposition (A = QR)");
    myMatrix.Transpose();
    myMatrix.Randomize();

    PrintInfo("The matrix A to decompose");
    std::cout << myMatrix << "\n\n";

    // QR Decomposition
    MatrixMN<double> Qthin, Rthin;
    myMatrix.GetQRDecomposition(Qthin, Rthin);

    PrintInfo("Q");
    std::cout << Qthin << "\n\n";
    PrintInfo("R");
    std::cout << Rthin << "\n\n";
    PrintInfo("Verification Q^t*Q = I)");
    MatrixMN<double> Q_transpose_mult_Q = Qthin.transpose() * Qthin;
    std::cout << Q_transpose_mult_Q << "\n\n";
    EXPECT_TRUE(Q_transpose_mult_Q.IsIdentity(1e-8));
    PrintInfo("Verification (A - Q*R)");
    std::cout << myMatrix - Qthin * Rthin << "\n\n";
    EXPECT_TRUE(myMatrix.IsEqual(Qthin * Rthin));

    PrintHeader("Full QR Decomposition (A = QR)");

    PrintInfo("The matrix A to decompose");
    std::cout << myMatrix << "\n\n";

    // QR Decomposition
    MatrixMN<double> Qfull, Rfull;
    myMatrix.GetFullQRDecomposition(Qfull, Rfull);

    PrintInfo("Q");
    std::cout << Qfull << "\n\n";
    PrintInfo("R");
    std::cout << Rfull << "\n\n";
    PrintInfo("Verification Q^t*Q = I)");
    MatrixMN<double> Q_transpose_mult_Q_02 = Qfull.transpose() * Qfull;
    std::cout << Q_transpose_mult_Q_02 << "\n\n";
    EXPECT_TRUE(Q_transpose_mult_Q_02.IsIdentity(1e-8));
    std::cout << Q_transpose_mult_Q_02 << "\n\n";
    PrintInfo("Verification (A - Q*R)");
    auto Rextended = MatrixMN<double>(6, 3);
    Rextended.SetNull();
    Rextended.block(0, 0, 3, 3) = Rfull;
    std::cout << myMatrix - Qfull * Rextended << "\n\n";
    EXPECT_TRUE(myMatrix.IsEqual(Qfull * Rextended));

    // ========================================================================
    //    LU Decomposition
    // ========================================================================

    PrintHeader("LU Decomposition (PA = LU) for A square");
    PrintInfo("The matrix A to decompose");
    matSquare.Randomize();
    std::cout << matSquare << "\n\n";

    // LU decomposition
    MatrixMN<double> P, L, U;
    matSquare.GetLUDecomposition(P, L, U);

    PrintInfo("P");
    std::cout << P << "\n\n";
    PrintInfo("L");
    std::cout << L << "\n\n";
    PrintInfo("U");
    std::cout << U << "\n\n";
    PrintInfo("Verification (P*A - L*U)");
    MatrixMN<double> PA = P * matSquare;
    std::cout << PA - L*U << "\n\n"; // Must be null !
    EXPECT_TRUE(PA.IsEqual(L*U));

    // ========================================================================
    //    Cholesky Decomposition
    // ========================================================================

    PrintHeader("Cholesky Decomposition A = L*L^T");
    PrintInfo("The matrix A to decompose");
    // Cholesky decomposition
    matSquare.Randomize();
    auto symMatrix = MakeSymmetricFromUpper(matSquare);
    symMatrix *= symMatrix;

    std::cout << symMatrix << "\n\n";

    MatrixMN<double> Lchol;
    symMatrix.GetCholeskyDecomposition(Lchol);

    PrintInfo("L");
    std::cout << Lchol << "\n\n";
//    std::cout << Transpose(Lchol) << "\n\n";
//    std::cout << myMatrix << "\n\n";
    PrintInfo("Verification (A - L * L^T)");
    std::cout << symMatrix - Lchol * Transpose(Lchol) << "\n\n";  // FIXME: ne donne pas le bon resultat !!!
    EXPECT_TRUE(symMatrix.IsEqual(Lchol * Transpose(Lchol)));
//    std::cout << myMatrix - L * Transpose(L) << "\n\n";

    auto myVector = MatrixMN<double>(6, 1);
    MatrixMN<double> LCholUpdate;
    matSquare.GetCholeskyUpdate(LCholUpdate,myVector,-1);


    PrintInfo("GetCholeskyUpdate:");
    std::cout << LCholUpdate << "\n\n";

    MatrixMN<double> LCholInplace(matSquare);
    LCholInplace.InPlaceCholeskyUpdate(myVector,-1);

    PrintInfo("GetCholeskyUpdate InPlace:");
    std::cout << LCholInplace << "\n\n";


    // ========================================================================
    //    SVD Decomposition
    // ========================================================================

    PrintHeader("SVD Decomposition A = U * Sdiag * (V*)");
    // SVD decomposition
    PrintInfo("The matrix to decompose");
    myMatrix.Randomize();
    std::cout << myMatrix << "\n\n";

    MatrixMN<double> Usvd, SingVal, Vsvd;
    myMatrix.GetSVDDecomposition(Usvd, SingVal, Vsvd);

    PrintInfo("U");
    std::cout << Usvd << "\n\n";
    PrintInfo("Singular values");
    std::cout << SingVal << "\n\n";
    PrintInfo("V");
    std::cout << Vsvd << "\n\n";

    PrintInfo("Verification (A - U * Sdiag * V*)");
    std::cout << myMatrix - Usvd * SingVal.asDiagonal() * Vsvd.adjoint() << "\n\n";
    EXPECT_TRUE(myMatrix.IsEqual(Usvd * SingVal.asDiagonal() * Vsvd.adjoint()));

    // ========================================================================
    //    Pseudo inverse computation
    // ========================================================================

    PrintHeader("Pseudo inverse of a rectangular matrix");
    // Pseudo inverse
    auto matRect = MatrixMN<double>(10, 5);
    matRect.Randomize();
    auto matRect_pinv = matRect.GetPseudoInverse();

    PrintInfo("A+");
    std::cout << matRect_pinv << "\n\n";

    PrintInfo("Using the function inversion");
    // Fonction Pinv
    auto pinvMat = Pinv(matRect);
    std::cout << pinvMat << "\n\n";

    PrintInfo("Verification (A+ * A)");
    MatrixMN<double> ApA = matRect_pinv * matRect;
    std::cout << ApA << "\n\n";

    EXPECT_TRUE(ApA.IsIdentity());

    // ========================================================================
    //    Linear system solvers
    // ========================================================================

    PrintHeader("Linear system solvers for Ax = b with A square");

    // Matrix A.
    auto matA = MatrixMN<double>(5, 5);
    matA.Randomize();
    PrintInfo("The matrix A");
    std::cout << matA << "\n\n";

    // The right-hand side vector.
    auto vectB = VectorN<double>(5);
    vectB.Randomize();
    PrintInfo("The right-hand side vector b");
    std::cout << vectB << "\n\n";

    // Solving using LU decomposition.
    auto solLU = matA.LUSolver(vectB);
    PrintInfo("The solution x from a LU decomposition");
    std::cout << solLU << "\n\n";

    // Solving using QR decomposition.
    auto solQR = matA.QRSolver(vectB);
    PrintInfo("The solution x from a QR decomposition");
    std::cout << solQR << "\n\n";

    PrintInfo("Verification (xLU - xQR = 0)");
    std::cout << solLU - solQR << "\n\n";
    EXPECT_TRUE(solLU.IsEqual(solQR));

    PrintInfo("Verification (A * xLU - b = 0)");
    std::cout << matA * solLU - vectB << "\n\n";
    VectorN<double> AxLu = matA * solLU;
    EXPECT_TRUE(AxLu.IsEqual(vectB));

    // The right-hand side matrix.
    auto matB = MatrixMN<double>(5, 2);
    matB.Randomize();
    PrintInfo("The right-hand side matrix B");
    std::cout << matB << "\n\n";

    // Solving using LU decomposition.
    auto solmat = matA.LUSolver(matB);
    PrintInfo("The solution x from with several right-hand side vectors");
    std::cout << solmat << "\n\n";

    PrintInfo("Verification with several right-hand side vectors (A * x - B = 0)");
    std::cout << matA * solmat - matB << "\n\n";
    MatrixMN<double> Ax = matA * solmat;
    EXPECT_TRUE(Ax.IsEqual(matB));

    // ========================================================================
    //    Linear least square system solver
    // ========================================================================

    PrintHeader("Linear least square system solver for Ax = b with A rectangular");

    // Matrix A.
    auto matALS = MatrixMN<double>(10, 5);
    matALS.Randomize();
    PrintInfo("The matrix A");
    std::cout << matALS << "\n\n";

    // The right-hand side vector.
    auto vectBLS = VectorN<double>(10);
    vectBLS.Randomize();
    PrintInfo("The right-hand side vector b");
    std::cout << vectBLS << "\n\n";

    // Solving using SVD decomposition.
    auto solLS = matALS.LeastSquareSolver(vectBLS);
    PrintInfo("The least square solution x from a SVD decomposition");
    std::cout << solLS << "\n\n";

    PrintInfo("Verification (x - pseudo_inv(A) * b = 0)");
    auto matALS_inv = matALS.GetPseudoInverse();
    std::cout << solLS - matALS_inv * vectBLS << "\n\n";
    VectorN<double> AxLS = matALS * solLS;
    EXPECT_TRUE(solLS.IsEqual(matALS_inv * vectBLS));
    // It is normal that A * x - b != 0.

    // ========================================================================
    //    Linear least square system solver under constraint
    // ========================================================================

    PrintHeader("Linear least square system solver for Ax = b with A rectangular subject to the constraint Cx = d.");

    // Matrix A.
    PrintInfo("The matrix A");
    std::cout << matALS << "\n\n";

    // The right-hand side vector.
    PrintInfo("The right-hand side vector b");
    std::cout << vectBLS << "\n\n";

    // Matrix C.
    auto matC = MatrixMN<double>(4, 5);
    matC.Randomize();
    PrintInfo("The matrix C");
    std::cout << matC << "\n\n";

    // The vector d.
    auto vectd = VectorN<double>(4);
    vectd.Randomize();
    PrintInfo("The vector d");
    std::cout << vectd << "\n\n";

    // Solving using SVD decomposition.
    auto solLSConstraint = matALS.LeastSquareSolverConstraint(vectBLS, matC, vectd);
    PrintInfo("The least square solution x subject to the constraint Cx = d");
    std::cout << solLSConstraint << "\n\n";

    PrintInfo("Verification (x - pseudo_inv(A) * b != 0)");
    std::cout << solLSConstraint - matALS_inv * vectBLS << std::endl;
    VectorN<double> AxLSConstraint = matALS * solLSConstraint;
    std::cout << "Because of the constraint, the equation x = inv(A) * b is not satisfied anymore. Some errors appear." << "\n\n";
    // It is normal that both A * x - b != 0 and x - pseudoinv(A) * b != 0.

    PrintInfo("Verification (Cx - d = 0)");
    std::cout << matC * solLSConstraint - vectd << "\n\n";
    VectorN<double> Cx = matC * solLSConstraint;
    EXPECT_TRUE(Cx.IsEqual(vectd));

    // ========================================================================
    //    Eigenvalues and eigenvectors.
    // ========================================================================

    PrintHeader("Eigenvalues and eigenvectors");

    // Matrix A.
    mathutils::MatrixMN<double> Avp(3, 3);
    Avp << 3,1,1,-8,-3,-4,6,3,4;
    PrintInfo("The matrix A");
    std::cout << Avp << "\n\n";

    // Computation of the eigenvalues.
    mathutils::MatrixMN<std::complex<double>> vp;
    vp = Avp.Eigenvalues();
    PrintInfo("The eigenvalues of A");
    std::cout << vp << "\n\n";
    PrintInfo("The eigenvalues must be (2, 1, 1)");

    int m_1 = 0, m_2 = 0;
    for (int i = 0; i < 3; i++) {
        double re = vp.at(i, 0).real();
        double im = vp.at(i, 0).imag();
        if (re < 1.5) {
            EXPECT_NEAR(re, 1.0, 1e-8);
            m_1 += 1;
        } else {
            EXPECT_NEAR(re, 2.0, 1e-8);
            m_2 += 1;
        }
        EXPECT_NEAR(im, 0.0, 1e-8);
    }
    EXPECT_EQ(m_1, 2);
    EXPECT_EQ(m_2, 1);
}
