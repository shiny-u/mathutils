//
// Created by frongere on 27/09/18.
//

#ifndef FRYDOM_MATRIX33_H
#define FRYDOM_MATRIX33_H

#include "EigenDense.h"
#include "iostream"

namespace mathutils {

    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================


    template <class Scalar_T>
    class Matrix33 : public Eigen::Matrix<Scalar_T, 3, 3> {

    public:

        using Scalar = Scalar_T;

        // =====================================================================
        //  Constructors
        // =====================================================================

        Matrix33() {};

//        Matrix33(unsigned int nbRows, unsigned int nbCols) : Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(nbRows, nbCols) {}

        // =====================================================================
        // Initialization methods
        // =====================================================================
        inline Scalar at(unsigned int irow, unsigned int icol) const;

        inline Scalar& at(unsigned int irow, unsigned int icol);

        void Randomize();

        void SetIdentity();

        void SetNull();

        void print(std::string name) const;

        // =====================================================================
        // Data Extraction methods
        // =====================================================================
        Eigen::Matrix<Scalar, 3, 1> GetColumn(unsigned int iCol) const;

        Eigen::Matrix<Scalar, 1, 3> GetRow(unsigned int iRow) const;

        Eigen::Matrix<Scalar, 3, 1> GetDiag() const;

        // =====================================================================
        // Matrix manipulation methods
        // =====================================================================
        void Transpose();

        // TODO: porter les methodes de resize, reshape...

        // =====================================================================
        // Matrix properties methods
        // =====================================================================

        unsigned int GetNbRows() const;

        unsigned int GetNbCols() const;

        bool IsPositiveSemiDefinite(const Scalar& epsilon=1e-7) const;

        bool IsSymmetric() const;

        bool IsIdentity() const;

        bool IsSquare() const;

        bool IsOrthogonal() const;

        bool IsEqual(const Matrix33<Scalar>& other, const Scalar& epsilon=1e-12) const;

        // =====================================================================
        // Various matrix decompositions
        // =====================================================================

        // As a matrix 3x3 is square, there is no need of GetFullQRDecomposition.
        void GetQRDecomposition(Matrix33<Scalar>& Q, Matrix33<Scalar>& R) const;

        void GetLUDecomposition(Matrix33<Scalar>& P, Matrix33<Scalar>& L, Matrix33<Scalar>& U) const;

        void GetCholeskyDecomposition(Matrix33<Scalar>& L) const;

        void GetSVDDecomposition(Matrix33<Scalar>& U, Matrix33<Scalar>& S, Matrix33<Scalar>& V) const;

        // =====================================================================
        // Various matrix inverse methods
        // =====================================================================
        void Inverse();

        Matrix33<Scalar> GetInverse() const;

        Matrix33<Scalar> GetPseudoInverse(Scalar tol=1e-6) const;

        // =====================================================================
        // Methods for Eigen inheritance usage
        // =====================================================================

        // This constructor allows to construct MatrixMN from Eigen expressions
        template <class OtherDerived>
        Matrix33(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, 3, 3>(other) {}

        // This method allows to assign Eigen expressions to MatrixMN
        template <class OtherDerived>
        Matrix33<Scalar>& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, 3, 3>::operator=(other);
            return *this;
        }

      // =====================================================================
      // Linear system solvers.
      // =====================================================================

      // A typename has to be used because the rhs is not necessary a Matrix33.
      template<typename T>
      T LUSolver(const T& rhs) const {
        return (this->fullPivLu().solve(rhs));
      }

      // A typename has to be used because the rhs is not necessary a Matrix33.
      template<typename T>
      T QRSolver(const T& rhs) const {
        return (this->fullPivHouseholderQr().solve(rhs));
      }

      // =====================================================================
      // Linear least square system solver.
      // =====================================================================

      // This method solved a least square problem min||Ax - b|| from a SVD decomposition (bidiagonal divide and
      // conquer SVD method).
      // A typename has to be used because the rhs is not necessary a Matrix66.
      template<typename T>
      T LeastSquareSolver(const T& b) const {
        return (this->bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b));
      }

      // No LeastSquareSolverContraint because it requires the structure MatrixMN for C and d, which is not available here.
      // Consequently, use a MatrixMN of size 3x3 for accessing this method.

      // =====================================================================
      // Eigenvalues and eigenvectors.
      // =====================================================================

      // This method computes the eigenvalues only.
      Matrix33<std::complex<double>> Eigenvalues() const;

    };

    // =================================================================================================================
    // =================================================================================================================
    //                                              FUNCTIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    Matrix33<Scalar> Transpose(const Matrix33<Scalar>& mat) {
        auto tMat = mat;
        tMat.Transpose();
        return tMat;
    }

    template <class Scalar>
    Matrix33<Scalar> Pinv(const Matrix33<Scalar>& mat, const Scalar tol=1e-6) {
        return mat.GetPseudoInverse(tol);
    }

    template <class Scalar>
    Matrix33<Scalar> MakeSymmetricFromUpper(const Matrix33<Scalar>& other) {
        Eigen::Matrix<Scalar, 3, 3> newMat;
        newMat = other.template selfadjointView<Eigen::Upper>();
        return newMat;
    }

    template <class Scalar>
    Matrix33<Scalar> MakeSymmetricFromLower(const Matrix33<Scalar>& other) {
        Eigen::Matrix<Scalar, 3, 3> newMat;
        newMat = other.template selfadjointView<Eigen::Lower>();
        return newMat;
    }


    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    inline Scalar Matrix33<Scalar>::at(unsigned int irow, unsigned int icol) const {
        return this->operator()(irow, icol);
    }

    template <class Scalar>
    inline Scalar& Matrix33<Scalar>::at(unsigned int irow, unsigned int icol) {
        return this->operator()(irow, icol);
    }

    template <class Scalar>
    unsigned int Matrix33<Scalar>::GetNbRows() const {
        return 3;
    }

    template <class Scalar>
    unsigned int Matrix33<Scalar>::GetNbCols() const {
        return 3;
    }

    template <class Scalar>
    void Matrix33<Scalar>::SetNull() {
        this->setZero();
    }

    template <class Scalar>
    void Matrix33<Scalar>::print(std::string name) const {
        std::cout << "\t" << name << ":\n";
        std::cout << std::endl << *this << std::endl;
    }

    template <class Scalar>
    void Matrix33<Scalar>::Randomize() {
        this->setRandom();
    }

    template <class Scalar>
    void Matrix33<Scalar>::SetIdentity() {
        this->setIdentity();
    }

    template <class Scalar>
    Eigen::Matrix<Scalar, 3, 1> Matrix33<Scalar>::GetColumn(unsigned int iCol) const {
        return this->col(iCol);
    }

    template <class Scalar>
    Eigen::Matrix<Scalar, 1, 3> Matrix33<Scalar>::GetRow(unsigned int iRow) const {
        return this->row(iRow);
    }

    template <class Scalar>
    Eigen::Matrix<Scalar, 3, 1> Matrix33<Scalar>::GetDiag() const {
        return this->diagonal();
    }

    template <class Scalar>
    void Matrix33<Scalar>::Transpose() {
        this->transposeInPlace();
    }

    template <class Scalar>
    void Matrix33<Scalar>::Inverse() {
        this->swap(this->inverse());
    }

    template <class Scalar>
    Matrix33<Scalar> Matrix33<Scalar>::GetInverse() const {
        auto newMat = *this;
        newMat.Inverse();
        return newMat;
    }

    template <class Scalar>
    bool Matrix33<Scalar>::IsPositiveSemiDefinite(const Scalar& epsilon) const {
        return (this->rows() == this->cols()) && ((this->eigenvalues().real().array() > -epsilon).all());
    }

    template <class Scalar>
    bool Matrix33<Scalar>::IsSymmetric() const {
        auto a = this->isApprox(this->adjoint());
        // FIXME: finir implementation :::
//        auto a = (this->adjoint() == this);
//        return this
//        if (this->adjoint() == this) {
//            return true;
//        }
        return true;
    }

    template <class Scalar>
    bool Matrix33<Scalar>::IsIdentity() const {
        return this->isIdentity();
    }

    template <class Scalar>
    bool Matrix33<Scalar>::IsSquare() const {
        return true;
    }

    template <class Scalar>
    bool Matrix33<Scalar>::IsOrthogonal() const {
        return this->isOrthogonal(*this);
    }

    template <class Scalar>
    bool Matrix33<Scalar>::IsEqual(const Matrix33<Scalar>& other, const Scalar& epsilon) const {
        return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    void Matrix33<Scalar>::GetQRDecomposition(Matrix33<Scalar>& Q, Matrix33<Scalar>& R) const {
        assert(this->rows() >= this->cols());

        auto QR = this->householderQr();

        // Q
        auto Qfull = QR.householderQ();
        auto thinQ = Matrix33<Scalar>::Identity(QR.rows(), QR.cols());
        Q = (Matrix33<Scalar>)(Qfull * thinQ);

        // R
        Eigen::Matrix<Scalar, 3, 3> Rfull;
        Rfull = QR.matrixQR().template triangularView<Eigen::Upper>();
        R = (Matrix33<Scalar>)(Rfull.block(0, 0, QR.cols(), QR.cols()));  // TODO: verifier que c'est bien cols() cols()...

    }

    template <class Scalar>
    void Matrix33<Scalar>::GetLUDecomposition(Matrix33<Scalar>& P, Matrix33<Scalar>& L, Matrix33<Scalar>& U) const {

        // PA = LU

        assert(IsSquare());

        auto LU = this->partialPivLu();

        // Permutation matrix P
        Eigen::Matrix<unsigned int, 3, 3> PP;
        PP = LU.permutationP();
        Eigen::Matrix<Scalar, 3, 3> PPScalar = PP.template cast<Scalar>();
        P = (Matrix33<Scalar>)(PPScalar);

        // L
        Eigen::Matrix<Scalar, 3, 3> Lpartial;
        Lpartial = LU.matrixLU().template triangularView<Eigen::StrictlyLower>();
        L = (Matrix33<Scalar>)(Lpartial + Eigen::Matrix<Scalar, 3, 3>::Identity(LU.rows(), LU.cols()));

        // U
        Eigen::Matrix<Scalar, 3, 3> Utmp;
        Utmp = LU.matrixLU().template triangularView<Eigen::Upper>();
        U = (Matrix33<Scalar>)(Utmp);

    }

    template <class Scalar>
    void Matrix33<Scalar>::GetCholeskyDecomposition(Matrix33<Scalar>& L) const {

        // Verifying that the matrix is square
        if (!IsSquare()) {
            throw std::runtime_error("In Cholesky decomposition, matrix must be square");
        }

        // Verifying that the matrix is symmetric
        if (!IsSymmetric()) {
            throw std::runtime_error("In Cholesky decomposition, matrix must be self adjoint (symmetric if real coefficients)");
        }

        auto CHOL = this->llt();

        if (CHOL.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("In Cholesky decomposition, matrix must be positive semi definite and this one is possibly not");
        }

        Eigen::Matrix<Scalar, 3, 3> Ltmp;
        Ltmp = CHOL.matrixL();
        L = (Matrix33<Scalar>)(Ltmp);
    }

    template <class Scalar>
    void Matrix33<Scalar>::GetSVDDecomposition(Matrix33<Scalar>& U, Matrix33<Scalar>& S, Matrix33<Scalar>& V) const {

        auto SVD = this->jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::Matrix<Scalar, 3, 3> Utmp, Stmp, Vtmp;

        // U
        Utmp = SVD.matrixU();
        U = (Matrix33<Scalar>)(Utmp);

        // S
        Stmp = SVD.singularValues();
        S = (Matrix33<Scalar>)(Stmp);

        // V
        Vtmp = SVD.matrixV();
        V = (Matrix33<Scalar>)(Vtmp);
    }

    template <class Scalar>
    Matrix33<Scalar> Matrix33<Scalar>::GetPseudoInverse(const Scalar tol) const {

        // SVD decomposition
        Matrix33<Scalar> U, S, V;
        GetSVDDecomposition(U, S, V);

        for (unsigned int i=0; i < S.rows(); ++i) {
            if (S(i) > tol) {
                S(i) = 1. / S(i);
            } else {
                S(i) = 0.;
            }
        }

        return V * S.asDiagonal() * U.adjoint();
    }


    template <class Scalar>
    Matrix33<std::complex<double>> Matrix33<Scalar>::Eigenvalues() const{

      // This method computes the eigenvalues only.

      // Verification.
      assert(this->GetNbRows() == this->GetNbCols());

      // Object to computing eigenvalues and eigenvectors.
      Eigen::EigenSolver<Eigen::MatrixXd> es(*this, false);

      return Matrix33<std::complex<double>>(es.eigenvalues());

    }

}  // end namesapce mathutils



#endif //FRYDOM_MATRIX33_H
