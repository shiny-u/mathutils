//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_VECTOR6D_H
#define MATHUTILS_VECTOR6D_H

#include "EigenDense.h"

#include "Unit.h"
#include "Matrix66.h"

namespace mathutils {
    // See the following link for Eigen::Matrix inheritance :
    // eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html

    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar_T=double>
    class Vector6d : public Eigen::Matrix<Scalar_T, 6, 1> {

    public:

        using Scalar = Scalar_T;

        Vector6d();

        Vector6d(Scalar x0, Scalar x1, Scalar x2, Scalar x3, Scalar x4, Scalar x5);

        // This constructor allows to construct Vector2d from Eigen expressions
        template <class OtherDerived>
        Vector6d(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, 6, 1>(other) {}

        // This method allows to assign Eigen expressions to Vector6d
        template <class OtherDerived>
        Vector6d& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, 6, 1>::operator=(other);
            return *this;
        }

        inline Scalar at(unsigned int index) const { return this->operator[](index); }

        inline Scalar& at(unsigned int index) { return this->operator[](index); }

        inline void SetNull();

        inline void Set(Scalar x0, Scalar x1, Scalar x2, Scalar x3, Scalar x4, Scalar x5);

        bool IsEqual(const Vector6d<Scalar>& other, const Scalar& epsilon=1e-12);

        inline Scalar infNorm() const;

        void print(std::string name) const;

    };

    // =================================================================================================================
    // Functions declarations
    // =================================================================================================================





    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    // =================================================================================================================
    // Vector6d methods implementations
    // =================================================================================================================

    template <class Scalar>
    Vector6d<Scalar>::Vector6d() : Eigen::Matrix<Scalar, 6, 1>() {
        SetNull();
    }

    template <class Scalar>
    Vector6d<Scalar>::Vector6d(Scalar x0, Scalar x1, Scalar x2, Scalar x3, Scalar x4, Scalar x5) {
        Set(x0, x1, x2, x3, x4, x5);
    }

    template <class Scalar>
    void Vector6d<Scalar>::SetNull() {
        Eigen::Matrix<Scalar, 6, 1>::setZero();
    }

    template <class Scalar>
    void Vector6d<Scalar>::Set(Scalar x0, Scalar x1, Scalar x2, Scalar x3, Scalar x4, Scalar x5) {
        this->operator[](0) = x0;
        this->operator[](1) = x1;
        this->operator[](2) = x2;
        this->operator[](3) = x3;
        this->operator[](4) = x4;
        this->operator[](5) = x5;
    }

    template <class Scalar>
    bool Vector6d<Scalar>::IsEqual(const Vector6d<Scalar>& other, const Scalar& epsilon) {
      return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    Scalar Vector6d<Scalar>::infNorm() const {
        return this->Eigen::Matrix<Scalar, 6, 1>::template lpNorm<Eigen::Infinity>();
    }

    template <class Scalar>
    void Vector6d<Scalar>::print(std::string name) const {
        std::cout << "\t" << name << ":\n";
        std::cout << std::endl << *this << std::endl;
    }

    // =================================================================================================================
    // Functions implementations
    // =================================================================================================================

    template <class Scalar>
    Matrix66<Scalar> outer(const Vector6d<Scalar>& v1, const Vector6d<Scalar>& v2) {
        return v1 * v2.transpose();  // TODO : verifier !!!
    }

    template <class Scalar>
    Matrix66<Scalar> outer(const Vector6d<Scalar>& v) {
        return outer(v, v);
    }

    template <class Scalar>
    Scalar inner(const Vector6d<Scalar>& v1, const Vector6d<Scalar>& v2) {
//        return v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) + v1(4)*v2(4) + v1(5)*v2(5);
        return v1.transpose() * v2;
    }

}  // end namespace mathutils


#endif //MATHUTILS_VECTOR6D_H
