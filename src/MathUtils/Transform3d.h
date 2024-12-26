//
// Created by frongere on 04/12/17.
//

#ifndef MATHUTILS_TRANSFORM3D_H
#define MATHUTILS_TRANSFORM3D_H

#include "Vector3d.h"

#include <Eigen/Geometry>

namespace mathutils {

    template <class Scalar>
    class Quaternion : Eigen::Quaternion<Scalar> {

    };

    template <class Scalar>
    class RotationMatrix : Eigen::Matrix<Scalar, 3, 3> {

    };


    template <class Scalar>
    class EulerAngleSequence : Vector3d<Scalar> {

    };


    template <class Scalar>
    class Rotation {

    private:
        Quaternion<Scalar> m_quaternion;
        RotationMatrix<Scalar> m_rotMat;

    public:
        Rotation();

        Rotation(const Quaternion<Scalar>& quat);

        void Set(const Quaternion<Scalar>& quat);
        void Set(const RotationMatrix<Scalar>& rotMat);
        void Set(const Vector3d<Scalar>& axis, const Scalar angle);
//        void Set(const Vector3d<Scalar>& eulerAngles, )

        Vector3d<Scalar> Apply(const Vector3d<Scalar>& vector3d) const;
        void Apply(Vector3d<Scalar>& vector3d) const;

        Vector3d<Scalar> InvApply(const Vector3d<Scalar>& vector3d) const;
        void InvApply(Vector3d<Scalar>& vector3d) const;


    };


    template <class Scalar>
    class Transform3d {

    private:
        Rotation<Scalar> m_rotation;
        Vector3d<Scalar> m_translation;



    };







}


#endif //MATHUTILS_TRANSFORM3D_H
