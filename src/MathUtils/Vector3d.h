//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_VECTOR3D_H
#define MATHUTILS_VECTOR3D_H

#include "EigenDense.h"

#include "Unit.h"
#include "Matrix33.h"
#include "Check.h"

namespace mathutils {
    // See the following link for Eigen::Matrix inheritance :
    // eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html

    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar_T=double>
    class Vector3d : public Eigen::Matrix<Scalar_T, 3, 1> {

    public:

        using Scalar = Scalar_T;

        Vector3d();

        Vector3d(Scalar x, Scalar y, Scalar z);

        // This constructor allows to construct Vector2d from Eigen expressions
        template <class OtherDerived>
        Vector3d(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, 3, 1>(other) {}

        // This method allows to assign Eigen expressions to Vector3d
        template <class OtherDerived>
        Vector3d& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, 3, 1>::operator=(other);
            return *this;
        }

        inline Scalar at(unsigned int index) const { return this->operator[](index); }

        inline Scalar& at(unsigned int index) { return this->operator[](index); }

        inline void SetNull();

        inline void Set(Scalar x, Scalar y, Scalar z);

        bool IsEqual(const Vector3d<Scalar>& other, const Scalar& epsilon=1e-12);

        inline Scalar infNorm() const;

        void Normalize();

        bool IsUnit() const;

        void print(std::string name) const;

        /// Return the angle between the y-axis and the vector projected in the YZ plane
        double GetProjectedAngleAroundX(ANGLE_UNIT unit) const {
            auto angle = atan2(this->at(2), this->at(1));
            if (unit == DEG) { angle *= RAD2DEG; }
            return angle;
        }

        /// Return the angle between the z-axis and the vector projected in the XZ plane
        double GetProjectedAngleAroundY(ANGLE_UNIT unit) const {
            auto angle = atan2(this->at(0), this->at(2));
            if (unit == DEG) { angle *= RAD2DEG; }
            return angle;
        }

        /// Return the angle between the x-axis and the vector projected in the XY plane
        double GetProjectedAngleAroundZ(ANGLE_UNIT unit) const {
            auto angle = atan2(this->at(1), this->at(0));
            if (unit == DEG) { angle *= RAD2DEG; }
            return angle;
        }
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
    // Vector3d methods implementations
    // =================================================================================================================

    template <class Scalar>
    Vector3d<Scalar>::Vector3d() : Eigen::Matrix<Scalar, 3, 1>() {
        SetNull();
    }

    template <class Scalar>
    Vector3d<Scalar>::Vector3d(Scalar x, Scalar y, Scalar z) {
        Set(x, y, z);
    }

    template <class Scalar>
    void Vector3d<Scalar>::SetNull() {
        Eigen::Matrix<Scalar, 3, 1>::setZero();
    }

    template <class Scalar>
    void Vector3d<Scalar>::Set(Scalar x, Scalar y, Scalar z) {
        this->operator[](0) = x;
        this->operator[](1) = y;
        this->operator[](2) = z;
    }

    template <class Scalar>
    bool Vector3d<Scalar>::IsEqual(const Vector3d<Scalar>& other, const Scalar& epsilon) {
      return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    Scalar Vector3d<Scalar>::infNorm() const {
        return this->Eigen::Matrix<Scalar, 3, 1>::template lpNorm<Eigen::Infinity>();
    }

    template <class Scalar>
    void Vector3d<Scalar>::Normalize() {
        this->normalize();
    }

    template <class Scalar>
    bool Vector3d<Scalar>::IsUnit() const {
        return IsClose<Scalar>(this->norm(), 1.);
    }

    template <class Scalar>
    void Vector3d<Scalar>::print(std::string name) const {
        std::cout << "\t" << name << ":\n";
        std::cout << std::endl << *this << std::endl;
    }

    // =================================================================================================================
    // Functions implementations
    // =================================================================================================================

    template <class Scalar>
    Matrix33<Scalar> outer(const Vector3d<Scalar>& v1, const Vector3d<Scalar>& v2) {
        return v1 * v2.transpose();
    }

    template <class Scalar>
    Matrix33<Scalar> outer(const Vector3d<Scalar>& v) {
        return outer(v, v);
    }

    template <class Scalar>
    double inner(const Vector3d<Scalar>& v1, const Vector3d<Scalar>& v2) {
        return v1.transpose() * v2;
    }

}  // end namespace mathutils


#endif //MATHUTILS_VECTOR3D_H
