//
// Created by frongere on 25/01/18.
//

#ifndef MATHUTILS_VECTORN_H
#define MATHUTILS_VECTORN_H

#include "EigenDense.h"
#include "Matrix.h"

namespace mathutils {

    // TODO: travailler sur les operateurs + et *
    // TODO: Rajouter methode Resize()


    template <class Scalar_T>
    class VectorN : public Eigen::Matrix<Scalar_T, Eigen::Dynamic, 1> {

    public:

        using Scalar = Scalar_T;

        VectorN() {}

        VectorN(unsigned int nbRows) : Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(nbRows) {}


        inline Scalar at(unsigned int n) const;

        inline Scalar& at(unsigned int n);

        void Randomize();

        void SetNull();

        void SetOnes();

        void Sort(bool ascending=true);

        void Reverse();

        bool IsEqual(const VectorN<Scalar>& other, const Scalar& epsilon=1e-12);

        MatrixMN<Scalar> GetMatrixSquare() const;

        MatrixMN<Scalar> GetDiagonalMatrix() const;

        inline Scalar infNorm() const;

        void print(std::string name) const;


        // =====================================================================
        // Methods for Eigen inheritance usage
        // =====================================================================

        // This constructor allows to construct VectorN from Eigen expressions
        template <class OtherDerived>
        VectorN(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(other) {}

        // This method allows to assign Eigen expressions to VectorN
        template <class OtherDerived>
        VectorN<Scalar>& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::operator=(other);
            return *this;
        }

    };


    // Functions

    template <class Scalar>
    MatrixMN<Scalar> GetMatrixSquare(const VectorN<Scalar> vector) {
        return vector * vector.transpose();
    }

    template <class Scalar>  // FIXME : Bug ?
    MatrixMN<Scalar> Transpose(const VectorN<Scalar> vector) {
        return vector.transpose();
    }

    template <class Scalar>
    MatrixMN<Scalar> GetDiagonalMatrix(const VectorN<Scalar> vector) {
      return vector.asDiagonal().toDenseMatrix();
    }


    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    inline Scalar VectorN<Scalar>::at(unsigned int n) const {
        return this->operator()(n);
    }

    template <class Scalar>
    inline Scalar& VectorN<Scalar>::at(unsigned int n) {
        return this->operator()(n);
    }

    template <class Scalar>
    void VectorN<Scalar>::Randomize() {
        this->setRandom();
    }

    template <class Scalar>
    void VectorN<Scalar>::SetNull() {
        this->setZero();
    }

    template <class Scalar>
    void VectorN<Scalar>::SetOnes() {
      this->setOnes();
    }

    template <class Scalar>
    void VectorN<Scalar>::Sort(bool ascending) {
        std::sort(this->data(), this->data()+this->size());

        if (!ascending) Reverse();

    }

    template <class Scalar>
    MatrixMN<Scalar> VectorN<Scalar>::GetMatrixSquare() const {
        return (*this) * ((*this).transpose());
    }

    template <class Scalar>
    MatrixMN<Scalar> VectorN<Scalar>::GetDiagonalMatrix() const {
      return (*this).asDiagonal().toDenseMatrix();
    }

    template <class Scalar>
    void VectorN<Scalar>::Reverse() {
        std::reverse(this->data(), this->data()+this->size());
    }

    template <class Scalar>
    bool VectorN<Scalar>::IsEqual(const VectorN<Scalar>& other, const Scalar& epsilon) {
        return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    Scalar VectorN<Scalar>::infNorm() const {
        return this->Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::template lpNorm<Eigen::Infinity>();
    }

    template <class Scalar>
    void VectorN<Scalar>::print(std::string name) const {
        std::cout << "\t" << name << ":\n";
        std::cout << std::endl << *this << std::endl;
    }

}  // End namespace mathutils



#endif //MATHUTILS_VECTORN_H
