//
// Created by frongere on 29/11/17.
//

#include "MathUtils/MathUtils.h"

#include <gtest/gtest.h>


using namespace mathutils;


template <class Scalar>
void print(std::string msg, Vector2d<Scalar> &vector) {
    std::cout << msg << std::endl
              << vector << std::endl << std::endl;
}


template <class Scalar>
void print(Vector2d<Scalar> &vector) {
    print("", vector);
}


TEST(Vector2d, OldMain) {
    // TODO Add assertions.

    Vector2d<double> a, b, c;

    print("By default, a a is initialized by 0 coefficients", a);

    a << 1, 2;
    print("Initializing a a with << operator", a);

    a[0] = 2;
    print("Setting a component with the [] operator", a);

    std::cout << "getting a component :" << std::endl << a[0] << std::endl << std::endl;

    a.SetNull();
    print("Setting coefficients to 0", a);


    a << 1, 2;
    b << 3, 4;
    print("Vector a", a);
    print("Vector b", b);

    c = a+b;
    print("a+b", c);

    c = a-b;
    print("a-b", c);

    a += b;
    print("a += b", a);

    a -= b;
    print("a -= b", a);

    a *= 2;
    print("a *= 2", a);

    a /= 2;
    print("a /= 2", a);

    print("Vector a", a);
    std::cout << "sum of coefficients of a\n" << a.sum() << "\n\n";
    std::cout << "product of coefficients of a\n" << a.prod() << "\n\n";
    std::cout << "mean of coefficients of a\n" << a.mean() << "\n\n";
    std::cout << "minimum of coefficients of a\n" << a.minCoeff() << "\n\n";
    std::cout << "maximum of coefficients of a\n" << a.maxCoeff() << "\n\n";

    std::cout << "squared norm of a\n" << a.squaredNorm() << "\n\n";
    std::cout << "Euclidian norm of a\n" << a.norm() << "\n\n";
    std::cout << "Infinity norm of a\n" << a.infNorm() << "\n\n";

    // New functions
    Vector2d<double> e;
    Vector2d<double> d(10, 20);
    print("Our constructor with coefficients initialization", d);

    d.at(0) = 0;
    print("Setting the first component with at() method", d);

    d.x() = 10;
    print("Setting the first component with the x() method", d);

    d.Set(20, 10);
    print("Using the Set() method", d);

    d.ProjectLocalOnNED(e, 10, DEG);
    print("ProjectLocalOnNED() copy version", e);

    Vector2d<double> w;
    e.ProjectNEDOnLocal(w, 10, DEG);
    print("ProjectNEDOnLocal() copy version", w);

    e.ProjectNEDOnLocal(10, DEG);
    print("ProjectNEDOnLocal() non copy version", e);

    auto north = ProjectOnNorthAxis<double>(d, 10, DEG);
    auto east = d.ProjectOnEastAxis(10, DEG);

    std::cout << "ProjectOnNorthAxis() function\n" << north << "\n\n";
    std::cout << "ProjectOnEastAxis() method\n" << east << "\n\n";

    c.Set(north, east);
    auto x = c.ProjectOnLocalXAxis(10, DEG);
    auto y = c.ProjectOnLocalYAxis(10, DEG);
    std::cout << "ProjectOnLocalXAxis()\n" << x << "\n\n";
    std::cout << "ProjectOnLocalYAxis()\n" << y << "\n\n";

    auto vel = GetNEDVelocityVector<double>(4, 0.5, 90, DEG, KNOT, KNOT);
    print("GetNEDVelocityVector<double>(4, 0.5, 90) -> 4 knot longi, 0.5 knot lat, at 90 deg", vel);

    Vector2d<double> velTransport;
    Vector2d<double> point(1, 0);

    vel.TransportAtPoint(velTransport, point, 2, DEGS);
    // FIXME: le point est la position relative du point B par rapport au point A dans le repere d'expression du vecteur...
    print("TransportAtPoint() copy version", velTransport);
}
