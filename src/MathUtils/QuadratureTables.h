//
// Created by pierre-yves on 14/01/20.
//

#ifndef MATHUTILS_QUADRATURETABLES_H
#define MATHUTILS_QUADRATURETABLES_H

#include <vector>

namespace mathutils {

  /**
  * Class for handling quadrature tables.
  */
  class QuadratureTables{

   public:

    /// Setter of the table.
    void Set_table(const std::vector<double>& alpha, const std::vector<double>& beta, const std::vector<double>& weight){
      this->m_alpha.push_back(alpha);
      this->m_beta.push_back(beta);
      this->m_weight.push_back(weight);
    }

    /// Getter of the alpha (area) coordinates of the Gauss points for a given order.
    std::vector<double> alpha(int order){
      return(m_alpha[order - 1]);
    }

    /// Getter of the beta (area) coordinates of the Gauss points for a given order.
    std::vector<double> beta(int order){
      return(m_beta[order - 1]);
    }

    /// Getter of the weights of the Gauss points for a given order.
    std::vector<double> weight(int order){
      return(m_weight[order - 1]);
    }

    /// This function returns the size of the tables.
    int get_size_tables(){
      return(m_weight.size());
    }

   private:

    /// Weights of the Gauss points.
    std::vector<std::vector<double>> m_weight;

    /// Alpha and beta (area) coordinates of the Gauss points.
    std::vector<std::vector<double>> m_alpha;
    std::vector<std::vector<double>> m_beta;

  };

}

#endif //MATHUTILS_QUADRATURETABLES_H
