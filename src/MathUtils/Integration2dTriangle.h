//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_INTEGRATION2DTRIANGLE_H
#define MATHUTILS_INTEGRATION2DTRIANGLE_H

#include "Integration2d.h"

namespace mathutils {

  /**
  * Class for computing 2d integrations over a triangle.
  */
  template<typename T>
  class Integration2dTriangle : public Integration2d<T> {

   public:

    /// Contructor of the class.
    Integration2dTriangle(Integrand<T>* F, const int &order)
        : Integration2d<T>(F, order) {

      // Quadrature tables.

      { // index = 0 - 1 Gauss point - order 1.
        std::vector<double> alpha = {1. / 3.};
        std::vector<double> beta = {1. / 3.};
        std::vector<double> weight = {1};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 1 - 3 Gauss points - order 2.
        std::vector<double> alpha = {0.166666666666670, 0.166666666666670, 0.666666666666670};
        std::vector<double> beta = {0.166666666666670, 0.666666666666670, 0.166666666666670};
        std::vector<double> weight = {0.333333333333330, 0.333333333333330, 0.333333333333330};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 2 - 4 Gauss point - order 3.
        std::vector<double> alpha = {0.333333333333330, 0.200000000000000, 0.200000000000000, 0.600000000000000};
        std::vector<double> beta = {0.333333333333330, 0.200000000000000, 0.600000000000000, 0.200000000000000};
        std::vector<double> weight = {-0.56250000000000, 0.520833333333330, 0.520833333333330, 0.520833333333330};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 3 - 6 Gauss points - order 4.
        std::vector<double> alpha = {0.445948490915970, 0.445948490915970, 0.108103018168070, 0.091576213509770,
                                     0.091576213509770, 0.816847572980460};
        std::vector<double> beta = {0.445948490915970, 0.108103018168070, 0.445948490915970, 0.091576213509770,
                                    0.816847572980460, 0.091576213509770};
        std::vector<double> weight = {0.223381589678010, 0.223381589678010, 0.223381589678010, 0.109951743655320,
                                      0.109951743655320, 0.109951743655320};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 4 - 7 Gauss points - order 5.
        std::vector<double> alpha = {0.333333333333330, 0.470142064105110, 0.470142064105110, 0.059715871789770,
                                     0.101286507323460, 0.101286507323460, 0.797426985353090};
        std::vector<double> beta = {0.333333333333330, 0.470142064105110, 0.059715871789770, 0.470142064105110,
                                    0.101286507323460, 0.797426985353090, 0.101286507323460};
        std::vector<double> weight = {0.225000000000000, 0.132394152788510, 0.132394152788510, 0.132394152788510,
                                      0.125939180544830, 0.125939180544830, 0.125939180544830};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 5 - 12 Gauss points - order 6.
        std::vector<double> alpha = {0.24928674517091, 0.24928674517091, 0.50142650965818, 0.06308901449150,
                                     0.06308901449150, 0.87382197101700, 0.31035245103378, 0.63650249912140,
                                     0.05314504984482, 0.63650249912140, 0.31035245103378, 0.05314504984482};
        std::vector<double> beta = {0.24928674517091, 0.50142650965818, 0.24928674517091, 0.06308901449150,
                                    0.87382197101700, 0.06308901449150, 0.63650249912140, 0.05314504984482,
                                    0.31035245103378, 0.31035245103378, 0.05314504984482, 0.63650249912140};
        std::vector<double> weight = {0.11678627572638, 0.11678627572638, 0.11678627572638, 0.05084490637021,
                                      0.05084490637021, 0.05084490637021, 0.08285107561837, 0.08285107561837,
                                      0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 6 - 13 Gauss points - order 7.
        std::vector<double> alpha = {0.33333333333333, 0.26034596607904, 0.26034596607904, 0.47930806784192,
                                     0.06513010290222, 0.06513010290222, 0.86973979419557, 0.31286549600487,
                                     0.63844418856981, 0.04869031542532, 0.63844418856981, 0.31286549600487,
                                     0.04869031542532};
        std::vector<double> beta = {0.33333333333333, 0.26034596607904, 0.47930806784192, 0.26034596607904,
                                    0.06513010290222, 0.86973979419557, 0.06513010290222, 0.63844418856981,
                                    0.04869031542532, 0.31286549600487, 0.31286549600487, 0.04869031542532,
                                    0.63844418856981};
        std::vector<double> weight = {-0.14957004446768, 0.17561525743321, 0.17561525743321, 0.17561525743321,
                                      0.05334723560884, 0.05334723560884, 0.05334723560884, 0.07711376089026,
                                      0.07711376089026, 0.07711376089026, 0.07711376089026, 0.07711376089026,
                                      0.07711376089026};
        this->m_tables->Set_table(alpha, beta, weight);
      }
      { // index = 7 - 16 Gauss points - order 8.
        std::vector<double> alpha = {0.33333333333333, 0.45929258829272, 0.45929258829272, 0.08141482341455,
                                     0.17056930775176, 0.17056930775176, 0.65886138449648, 0.05054722831703,
                                     0.05054722831703, 0.89890554336594, 0.26311282963464, 0.72849239295540,
                                     0.00839477740996, 0.72849239295540, 0.26311282963464, 0.00839477740996};
        std::vector<double> beta = {0.33333333333333, 0.45929258829272, 0.08141482341455, 0.45929258829272,
                                    0.17056930775176, 0.65886138449648, 0.17056930775176, 0.05054722831703,
                                    0.89890554336594, 0.05054722831703, 0.72849239295540, 0.00839477740996,
                                    0.26311282963464, 0.26311282963464, 0.00839477740996, 0.72849239295540};
        std::vector<double> weight = {0.14431560767779, 0.09509163426728, 0.09509163426728, 0.09509163426728,
                                      0.10321737053472, 0.10321737053472, 0.10321737053472, 0.03245849762320,
                                      0.03245849762320, 0.03245849762320, 0.02723031417443, 0.02723031417443,
                                      0.02723031417443, 0.02723031417443, 0.02723031417443, 0.02723031417443};
        this->m_tables->Set_table(alpha, beta, weight);
      }

      // Warning.
      if (this->m_order > this->m_tables->get_size_tables()) {
        std::cout << "\nToo high order of quadrature for triangle. Use lower order." << std::endl;
        exit(0);
      }
      if (this->m_order < 1) {
        std::cout << "\nToo low order of quadrature for triangle. Use higher order." << std::endl;
        exit(0);
      }

    }

    /// This function computes the 2d integration over a triangle.
    T Compute(const Vector3d<double> &vertex_1, const Vector3d<double> &vertex_2,
                 const Vector3d<double> &vertex_3) const override {

      // Face area.
      Vector3d<double> e1 = vertex_1 - vertex_3;
      Vector3d<double> e2 = vertex_2 - vertex_3;
      double area = 0.5 * ((e1.cross(e2))).norm();

      // Computation of the 2d integration over a standard triangle.
      T integral = ComputeOverStandardTriangle(vertex_1, vertex_2, vertex_3);

      // 2d integration over the real triangle.
      integral *= 2*area;

      return(integral);

    }

   private:

    /// This function computes the 2d integration over a standard triangle.
    T ComputeOverStandardTriangle(const Vector3d<double> &vertex_1, const Vector3d<double> &vertex_2,
                                     const Vector3d<double> &vertex_3) const {

      // Quadrature coefficients.
      std::vector<double> alpha = this->m_tables->alpha(this->m_order);
      std::vector<double> beta = this->m_tables->beta(this->m_order);
      std::vector<double> weight = this->m_tables->weight(this->m_order);

      // Computation of the 2d integration.
      T result;
      result *= 0; // This works for scalars, vectors and matrices.
      for (unsigned int i = 0; i < weight.size(); i++) {

        // Position of the Gauss points.
        Vector3d<double> x = (1 - alpha.at(i) - beta.at(i)) * vertex_1 + alpha.at(i) * vertex_2 + beta.at(i) * vertex_3;

        // When T = VectorN, the size of result must be specified otherwise the += tmp_val does not work.
        // This is done by Evaluate for i = 0.
        // There is always i = 0 because m_tables has a minimum order of 1.
        if(i == 0){
          result = weight.at(i) * this->m_integrand->Evaluate(x);
        }
        else {
          T tmp_val = this->m_integrand->Evaluate(x);
          tmp_val *= weight.at(i);
          result += tmp_val;
        }

      }
      result *= 0.5; // The 1/2 coefficient matches the area of the standard triangle.

      return (result);

    }

  };

}

#endif //MATHUTILS_INTEGRATION2DTRIANGLE_H
