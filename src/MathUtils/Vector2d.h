//
// Created by frongere on 15/11/17.
//

#ifndef MATHUTILS_VECTOR2D_H
#define MATHUTILS_VECTOR2D_H

#include "EigenDense.h"

#include "Unit.h"
#include "Angles.h"

namespace mathutils {
    // See the following link for Eigen::Matrix inheritance :
    // eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html

    // =================================================================================================================
    // =================================================================================================================
    //                                              DECLARATIONS
    // =================================================================================================================
    // =================================================================================================================

    // =================================================================================================================
    // Class Vector2d declaration
    // =================================================================================================================

    template <class Scalar_T=double>
    class Vector2d : public Eigen::Matrix<Scalar_T, 2, 1> {

    public:

        using Scalar = Scalar_T;

        Vector2d();

        Vector2d(Scalar x, Scalar y);

        inline Scalar at(unsigned int index) const { return this->operator[](index); }

        inline Scalar& at(unsigned int index) { return this->operator[](index); }

        inline void SetNull();

        inline void Set(Scalar x, Scalar y);

        inline void SetNEDFromLocal(Scalar x, Scalar y, Scalar angle, ANGLE_UNIT unit=RAD);

        bool IsEqual(const Vector2d<Scalar>& other, const Scalar& epsilon=1e-12);

        inline Scalar infNorm() const;

        inline void ProjectLocalOnNED(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit=RAD) const;
        inline void ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit=RAD); /// Non const version (no copy)

        inline void ProjectNEDOnLocal(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit=RAD) const;
        inline void ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit=RAD); /// Non const version (no copy)

        inline Scalar ProjectOnNorthAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar ProjectOnEastAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar ProjectOnLocalXAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline Scalar ProjectOnLocalYAxis(Scalar angle, ANGLE_UNIT unit=RAD) const;

        inline void TransportAtPoint(Vector2d<Scalar>& out, const Vector2d<Scalar>& point,
                                     Scalar w, FREQUENCY_UNIT funit=RADS) const;
        inline void TransportAtPoint(const Vector2d<Scalar>& point,
                                     Scalar w, FREQUENCY_UNIT funit=RADS); /// Non const version (no copy)
        void print(std::string name) const;

        // ====================================
        // Methods for Eigen inheritance usage
        // ====================================

        // This constructor allows to construct Vector2d from Eigen expressions
        template <class OtherDerived>
        Vector2d(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<Scalar, 2, 1>(other) {}

        // This method allows to assign Eigen expressions to Vector2d
        template <class OtherDerived>
        Vector2d& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            this->Eigen::Matrix<Scalar, 2, 1>::operator=(other);
            return *this;
        }

    };

    // =================================================================================================================
    // Functions declarations
    // =================================================================================================================
    template <class Scalar=double>
    Vector2d<Scalar> GetNEDVelocityVector(Scalar plongitudinalSpeed, Scalar plateralSpeed, Scalar angle,
                                          ANGLE_UNIT angleUnit=RAD,
                                          SPEED_UNIT inputSpeedUnit=KNOT, SPEED_UNIT outputSpeedUnit=MS);

    template <class Scalar>
    Vector2d<Scalar> ProjectLocalOnNED(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Vector2d<Scalar> ProjectNEDOnLocal(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnNorthAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnEastAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnLocalXAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Scalar ProjectOnLocalYAxis(const Vector2d<Scalar>& vector, Scalar angle, ANGLE_UNIT unit=RAD);

    template <class Scalar>
    Vector2d<Scalar> TransportAtPoint(const Vector2d<Scalar> vector, const Vector2d<Scalar> point,
                                      Scalar w, FREQUENCY_UNIT funit=RADS);


    // =================================================================================================================
    // =================================================================================================================
    //                                              IMPLEMENTATIONS
    // =================================================================================================================
    // =================================================================================================================

    // =================================================================================================================
    // Vector2d methods implementations
    // =================================================================================================================

    template <class Scalar>
    Vector2d<Scalar>::Vector2d() : Eigen::Matrix<Scalar, 2, 1>() {
        SetNull();
    }

    template <class Scalar>
    Vector2d<Scalar>::Vector2d(Scalar x, Scalar y) {
        Set(x, y);
    }

    template <class Scalar>
    void Vector2d<Scalar>::Set(Scalar x, Scalar y) {
        this->operator[](0) = x;
        this->operator[](1) = y;
    }

    template <class Scalar>
    void Vector2d<Scalar>::SetNEDFromLocal(Scalar x, Scalar y, Scalar angle, ANGLE_UNIT unit){
        Set(x,y);
        ProjectLocalOnNED(angle, unit);
    }


    template <class Scalar>
    void Vector2d<Scalar>::SetNull() {
        Eigen::Matrix<Scalar, 2, 1>::setZero();
    }

    template <class Scalar>
    bool Vector2d<Scalar>::IsEqual(const Vector2d<Scalar>& other, const Scalar& epsilon) {
      return this->isApprox(other, epsilon);
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectLocalOnNED(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit) const {
        out = *this;
        out.ProjectLocalOnNED(angle, unit);
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectLocalOnNED(Scalar angle, ANGLE_UNIT unit) {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        auto x = c * at(0) - s * at(1);
        at(1) = s * at(0) + c * at(1);
        at(0) = x;
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectNEDOnLocal(Vector2d<Scalar>& out, Scalar angle, ANGLE_UNIT unit) const {
        out = *this;
        out.ProjectNEDOnLocal(angle, unit);
    }

    template <class Scalar>
    void Vector2d<Scalar>::ProjectNEDOnLocal(Scalar angle, ANGLE_UNIT unit) {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        auto x =  c * at(0) + s * at(1);
        at(1) = -s * at(0) + c * at(1);
        at(0) = x;
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnNorthAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return c * at(0) - s * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnEastAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return s * at(0) + c * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnLocalXAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return c * at(0) + s * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::ProjectOnLocalYAxis(Scalar angle, ANGLE_UNIT unit) const {
        Scalar c, s;
        _GetTrigo(angle, c, s, unit);
        return -s * at(0) + c * at(1);
    }

    template <class Scalar>
    Scalar Vector2d<Scalar>::infNorm() const {
        return this->Eigen::Matrix<Scalar, 2, 1>::maxCoeff();
    }

    template <class Scalar>
    inline void Vector2d<Scalar>::TransportAtPoint(Vector2d<Scalar>& out, const Vector2d<Scalar>& point,
                                                   Scalar pw, FREQUENCY_UNIT funit) const {
        out = *this;
        out.TransportAtPoint(point, pw, funit);
    }

    template <class Scalar>
    inline void Vector2d<Scalar>::TransportAtPoint(const Vector2d<Scalar>& point,
                                                   Scalar pw, FREQUENCY_UNIT funit) {
        auto w = pw;
        if (funit != RADS) {
            w = convert_frequency(w, funit, RADS);
        }

        at(0) = at(0) - point[1] * w;
        at(1) = at(1) + point[0] * w;
    }

    template <class Scalar>
    void Vector2d<Scalar>::print(std::string name) const {
        std::cout << "\t" << name << ":\n";
        std::cout << std::endl << *this << std::endl;
    }

    // =================================================================================================================
    // Functions implementations
    // =================================================================================================================

    template<class Scalar>
    Vector2d<Scalar> ProjectLocalOnNED(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        Vector2d<Scalar> out(vector);
        out.ProjectLocalOnNED(angle, unit);
        return out;
    }

    template<class Scalar>
    Scalar ProjectOnNorthAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnNorthAxis(angle, unit);
    }

    template<class Scalar>
    Vector2d<Scalar> ProjectNEDOnLocal(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        Vector2d<Scalar> out = vector;
        out.ProjectNEDOnLocal(angle, unit);
        return out;
    }

    template<class Scalar>
    Scalar ProjectOnEastAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnEastAxis(angle, unit);
    }

    template<class Scalar>
    Scalar ProjectOnLocalXAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnLocalXAxis(angle, unit);
    }

    template<class Scalar>
    Scalar ProjectOnLocalYAxis(const Vector2d<Scalar> &vector, Scalar angle, ANGLE_UNIT unit) {
        return vector.ProjectOnLocalYAxis(angle, unit);
    }

    template<class Scalar>
    Vector2d<Scalar>
    GetNEDVelocityVector(Scalar plongitudinalSpeed, Scalar plateralSpeed, Scalar angle, ANGLE_UNIT angleUnit,
                         SPEED_UNIT inputSpeedUnit, SPEED_UNIT outputSpeedUnit) {

        Scalar longitudinalSpeed = plongitudinalSpeed;
        Scalar lateralSpeed      = plateralSpeed;

        if (inputSpeedUnit != outputSpeedUnit) {
            longitudinalSpeed = convert_velocity_unit(longitudinalSpeed, inputSpeedUnit, outputSpeedUnit);
            lateralSpeed      = convert_velocity_unit(lateralSpeed, inputSpeedUnit, outputSpeedUnit);
        }

        auto localVect = Vector2d<Scalar>(longitudinalSpeed, lateralSpeed);
        return ProjectLocalOnNED(localVect, angle, angleUnit);
//        return outVect.ProjectLocalOnNED(angle, angleUnit);  // TODO: verifier qu'on appelle bien la version inline
    }

    template <class Scalar>
    Vector2d<Scalar> TransportAtPoint(const Vector2d<Scalar> vector, const Vector2d<Scalar> point,
                                      Scalar w, FREQUENCY_UNIT funit) {
        Vector2d<Scalar> out(vector);
        out.TransportAtPoint(point, w, funit);
        return out;
    }


}  // end namespace mathutils

#endif //MATHUTILS_VECTOR2D_H
