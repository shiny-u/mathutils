//
// Created by frongere on 14/12/17.
//

#ifndef MATHUTILS_MATRIX_H
#define MATHUTILS_MATRIX_H

#include "EigenDense.h"
#include "iostream"

namespace mathutils {
    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================


    template <class Scalar_T>
    class MatrixMN : public Eigen::Matrix<Scalar_T, Eigen::Dynamic, Eigen::Dynamic> {

    public:

        using Scalar = Scalar_T;

        // =====================================================================
        //  Constructors
        // =====================================================================

        MatrixMN() {};

        MatrixMN(unsigned int nbRows, unsigned int nbCols) : Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(nbRows, nbCols) {}

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
        MatrixMN<Scalar> GetColumn(unsigned int iCol) const;

//        MatrixMN<Scalar>& GetColumn(unsigned int iCol); // TODO: voir pour des methodes non const pour changer les cols/rows

        MatrixMN<Scalar> GetRow(unsigned int iRow) const;

//        MatrixMN<Scalar>& GetRow(unsigned int iRow);

        MatrixMN<Scalar> GetDiag() const;

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

        bool IsIdentity(const Scalar& epsilon) const;

        bool IsSquare() const;

        bool IsOrthogonal() const;

        bool IsEqual(const MatrixMN<Scalar>& other, const Scalar& epsilon=1e-12) const;

        // =====================================================================
        // Various matrix decompositions
        // =====================================================================

        // This methods applies a thin QR decomposition, for more information cf the technical notes.
        // If your matrix is square, GetQRDecomposition and GetFullQRDecomposition will give the same result with
        // the same performances.
        void GetQRDecomposition(MatrixMN<Scalar>& Q, MatrixMN<Scalar>& R) const;

        // This methods applies a full QR decomposition, for more information cf the technical notes.
        // If your matrix is square, GetQRDecomposition and GetFullQRDecomposition will give the same result with
        // the same performances.
        void GetFullQRDecomposition(MatrixMN<Scalar>& Q, MatrixMN<Scalar>& R) const;

        void GetLUDecomposition(MatrixMN<Scalar>& P, MatrixMN<Scalar>& L, MatrixMN<Scalar>& U) const;

        void GetCholeskyDecomposition(MatrixMN<Scalar>& L) const;

        void GetCholeskyUpdate(MatrixMN<Scalar>& L, MatrixMN<Scalar>& v, const Scalar & sigma) const;

        void InPlaceCholeskyUpdate(MatrixMN<Scalar>& v, const Scalar & sigma);

        void GetSVDDecomposition(MatrixMN<Scalar>& U, MatrixMN<Scalar>& S, MatrixMN<Scalar>& V) const;

        // =====================================================================
        // Various matrix inverse methods
        // =====================================================================
        void Inverse();

        MatrixMN<Scalar> GetInverse() const;

        MatrixMN<Scalar> GetPseudoInverse(Scalar tol=1e-6) const;

        // =====================================================================
        // Methods for Eigen inheritance usage
        // =====================================================================

        // This constructor allows to construct MatrixMN from Eigen expressions
        template <class OtherDerived>
        MatrixMN(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(other) {}

        // This method allows to assign Eigen expressions to MatrixMN
        template <class OtherDerived>
        MatrixMN<Scalar>& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::operator=(other);
            return *this;
        }

      // =====================================================================
      // Linear system solvers.
      // =====================================================================

      MatrixMN<Scalar> LUSolver(const MatrixMN<Scalar>& rhs) const;

      MatrixMN<Scalar> QRSolver(const MatrixMN<Scalar>& rhs) const;

      // =====================================================================
      // Linear least square system solvers.
      // =====================================================================

      // This method solved a least square problem min||Ax - b|| from a SVD decomposition (bidiagonal divide and
      // conquer SVD method).
      MatrixMN<Scalar> LeastSquareSolver(const MatrixMN<Scalar>& b) const;

      // This method solves a least square problem min||Ax - b||^2 subject to the equality constraint Cx = d.
      MatrixMN<Scalar> LeastSquareSolverConstraint(const MatrixMN<Scalar>& b, const MatrixMN<Scalar>& C
          , const MatrixMN<Scalar>& d) const;

      // =====================================================================
      // Eigenvalues and eigenvectors.
      // =====================================================================

      // This method computes the eigenvalues only.
      MatrixMN<std::complex<double>> Eigenvalues() const;

    };

    // =================================================================================================================
    // =================================================================================================================
    //                                              FUNCTIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    MatrixMN<Scalar> Transpose(const MatrixMN<Scalar>& mat) {
        auto tMat = mat;
        tMat.Transpose();
        return tMat;
    }

    template <class Scalar>
    MatrixMN<Scalar> Pinv(const MatrixMN<Scalar>& mat, const Scalar tol=1e-6) {
        return mat.GetPseudoInverse(tol);
    }

    template <class Scalar>
    MatrixMN<Scalar> MakeSymmetricFromUpper(const MatrixMN<Scalar>& other) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> newMat;
        newMat = other.template selfadjointView<Eigen::Upper>();
        return newMat;
    }

    template <class Scalar>
    MatrixMN<Scalar> MakeSymmetricFromLower(const MatrixMN<Scalar>& other) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> newMat;
        newMat = other.template selfadjointView<Eigen::Lower>();
        return newMat;
    }


    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    inline Scalar MatrixMN<Scalar>::at(unsigned int irow, unsigned int icol) const {
        // return this->operator()(irow * (this->cols()-1) + icol);
        return this->operator()(irow, icol);
    }

    template <class Scalar>
    inline Scalar& MatrixMN<Scalar>::at(unsigned int irow, unsigned int icol) {
        // return this->operator()(irow * (this->cols()-1) + icol);
        return this->operator()(irow, icol);
    }

    template <class Scalar>
    unsigned int MatrixMN<Scalar>::GetNbRows() const {
        return (unsigned int)this->rows();
    }

    template <class Scalar>
    unsigned int MatrixMN<Scalar>::GetNbCols() const {
        return (unsigned int)this->cols();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::SetNull() {
        this->setZero();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::print(std::string name) const {
        std::cout << "\t" << name << ":\n";
        std::cout << std::endl << *this << std::endl;
    }

    template <class Scalar>
    void MatrixMN<Scalar>::Randomize() {
        this->setRandom();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::SetIdentity() {
        this->setIdentity();
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetColumn(unsigned int iCol) const {
        return this->col(iCol);
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetRow(unsigned int iRow) const {
        return this->row(iRow);
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetDiag() const {
        return this->diagonal();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::Transpose() {
        this->transposeInPlace();
    }

    template <class Scalar>
    void MatrixMN<Scalar>::Inverse() {
        this->swap(this->inverse());
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetInverse() const {
        auto newMat = *this;
        newMat.Inverse();
        return newMat;
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsPositiveSemiDefinite(const Scalar& epsilon) const {
        return (this->rows() == this->cols()) && ((this->eigenvalues().real().array() > -epsilon).all());
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsSymmetric() const {
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
    bool MatrixMN<Scalar>::IsIdentity() const {
        return this->isIdentity();
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsIdentity(const Scalar& epsilon) const {
        return this->isIdentity(epsilon);
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsSquare() const {
        return (this->rows() == this->cols());
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsOrthogonal() const {
        return this->isOrthogonal(*this);
    }

    template <class Scalar>
    bool MatrixMN<Scalar>::IsEqual(const MatrixMN<Scalar>& other, const Scalar& epsilon) const {
        return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    void MatrixMN<Scalar>::GetQRDecomposition(MatrixMN<Scalar>& Q, MatrixMN<Scalar>& R) const {

        assert(this->rows() >= this->cols());

        auto QR = this->householderQr();

        // Qthin
        auto Qfull = QR.householderQ();
        auto thinQ = MatrixMN<Scalar>::Identity(QR.rows(), QR.cols());
        Q = (MatrixMN<Scalar>)(Qfull * thinQ);

        // R
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Rfull;
        Rfull = QR.matrixQR().template triangularView<Eigen::Upper>();
        R = (MatrixMN<Scalar>)(Rfull.block(0, 0, QR.cols(), QR.cols()));  // TODO: verifier que c'est bien cols() cols()...

    }

  template <class Scalar>
  void MatrixMN<Scalar>::GetFullQRDecomposition(MatrixMN<Scalar>& Q, MatrixMN<Scalar>& R) const {

      assert(this->rows() >= this->cols());

      auto QR = this->householderQr();

      // Qfull
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Qfull = QR.householderQ();
      Q = (MatrixMN<Scalar>)(Qfull);

      // R
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Rfull;
      Rfull = QR.matrixQR().template triangularView<Eigen::Upper>();
      R = (MatrixMN<Scalar>)(Rfull.block(0, 0, QR.cols(), QR.cols()));  // TODO: verifier que c'est bien cols() cols()...

  }

    template <class Scalar>
    void MatrixMN<Scalar>::GetLUDecomposition(MatrixMN<Scalar>& P, MatrixMN<Scalar>& L, MatrixMN<Scalar>& U) const {

        // PA = LU

        assert(IsSquare());

        auto LU = this->partialPivLu();

        // Permutation matrix P
        Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> PP;
        PP = LU.permutationP();
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> PPScalar = PP.template cast<Scalar>();
        P = (MatrixMN<Scalar>)(PPScalar);

        // L
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Lpartial;
        Lpartial = LU.matrixLU().template triangularView<Eigen::StrictlyLower>();
        L = (MatrixMN<Scalar>)(Lpartial + Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(LU.rows(), LU.cols()));

        // U
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Utmp;
        Utmp = LU.matrixLU().template triangularView<Eigen::Upper>();
        U = (MatrixMN<Scalar>)(Utmp);

    }

    template <class Scalar>
    void MatrixMN<Scalar>::GetCholeskyDecomposition(MatrixMN<Scalar>& L) const {

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

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Ltmp;
        Ltmp = CHOL.matrixL();
        L = (MatrixMN<Scalar>)(Ltmp);
    }


    template <class Scalar>
    void MatrixMN<Scalar>::GetCholeskyUpdate(MatrixMN<Scalar>& L, MatrixMN<Scalar>& v,const Scalar & sigma) const
    {
        auto CHOL = this->llt();
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Ltmp = CHOL.rankUpdate(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(v),sigma).matrixL();
        L = (MatrixMN<Scalar>)(Ltmp);
    }


    template <class Scalar>
    void MatrixMN<Scalar>::InPlaceCholeskyUpdate(MatrixMN<Scalar>& v, const Scalar & sigma) {
        ///TODO : improve by using dedicated inplace methods from Eigen
        auto CHOL = this->llt();
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Ltmp = CHOL.rankUpdate(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(v),sigma).matrixL();
        MatrixMN<Scalar> out = (MatrixMN<Scalar>)(Ltmp);
        this->swap(out);
    }


    template <class Scalar>
    void MatrixMN<Scalar>::GetSVDDecomposition(MatrixMN<Scalar>& U, MatrixMN<Scalar>& S, MatrixMN<Scalar>& V) const {

        auto SVD = this->jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Utmp, Stmp, Vtmp;

        // U
        Utmp = SVD.matrixU();
        U = (MatrixMN<Scalar>)(Utmp);

        // S
        Stmp = SVD.singularValues();
        S = (MatrixMN<Scalar>)(Stmp);

        // V
        Vtmp = SVD.matrixV();
        V = (MatrixMN<Scalar>)(Vtmp);
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::GetPseudoInverse(const Scalar tol) const {

        // SVD decomposition
        MatrixMN<Scalar> U, S, V;
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
    MatrixMN<Scalar> MatrixMN<Scalar>::LUSolver(const MatrixMN<Scalar>& rhs) const {
      return (this->fullPivLu().solve(rhs));
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::QRSolver(const MatrixMN<Scalar>& rhs) const {
      return (this->fullPivHouseholderQr().solve(rhs));
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::LeastSquareSolver(const MatrixMN<Scalar>& b) const {
      return (this->bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b));
    }

    template <class Scalar>
    MatrixMN<Scalar> MatrixMN<Scalar>::LeastSquareSolverConstraint(const MatrixMN<Scalar>& b, const MatrixMN<Scalar>& C
      , const MatrixMN<Scalar>& d) const {

        // This method solves a least square problem min||Ax - b||^2 subject to the equality constraint Cx = d.

        // Number of rows in A.
        int m = this->GetNbRows();

        // Number of columns in A.
        int n = this->GetNbCols();

        // Size of d.
        int p = d.GetNbRows();

        // Verifications of the sizes.
        assert(n > p);
        assert(m > n);
        assert(b.GetNbRows() == m);
        assert(C.GetNbCols() == n);
        assert(C.GetNbRows() == p);

        // Only a single vector for the right-hand sides.
        assert(b.GetNbCols() == 1);
        assert(d.GetNbCols() == 1);

        // QR factorisation of transpose(C).
        MatrixMN<Scalar> Q, R;
        MatrixMN<Scalar> Ct = C;
        Ct.Transpose();
        Ct.GetFullQRDecomposition(Q, R);

        // Solving transpose(R) * y = d.
        R.Transpose(); // From now on, in the variable R there is the transpose of R.
        auto y = R.LUSolver(d);

        // Creation of A1 and A2.
        MatrixMN<Scalar> AQ = *this * Q;
        MatrixMN<Scalar> A1 = AQ.block(0, 0, m, p);
        MatrixMN<Scalar> A2 = AQ.block(0, p, m, n - p);

        // LS problem wrt z.
        auto z = A2.LeastSquareSolver(b - A1 * y); // R has been transposed.

        // Solution.
        MatrixMN<Scalar> vect_y_z = MatrixMN<Scalar>(n, 1);
        vect_y_z.block(0, 0, p, 1) = y;
        vect_y_z.block(p, 0, n - p, 1) = z;
        auto x = Q * vect_y_z;

        return x;

    };

    template <class Scalar>
    MatrixMN<std::complex<double>> MatrixMN<Scalar>::Eigenvalues() const{

        // This method computes the eigenvalues only.

        // Verification.
        assert(this->GetNbRows() == this->GetNbCols());

        // Object to computing eigenvalues and eigenvectors.
        Eigen::EigenSolver<Eigen::MatrixXd> es(*this, false);

        return MatrixMN<std::complex<double>>(es.eigenvalues());

    }


}  // end namespace mathutils

#endif //MATHUTILS_MATRIX_H

