//
// Created by frongere on 05/01/2021.
//

#ifndef MATHUTILS_MULTILINEARINTERPOLATOR_H
#define MATHUTILS_MULTILINEARINTERPOLATOR_H

#include <vector>
#include <array>
#include <memory>
#include <utility>
#include "thirdparty/linterp.h"

namespace mathutils {

// TODO: faire disparaitre cette fonction, cf la methode _BuildInterpolator (copier coler doc linterp)...
  template<class IterT>
  std::pair<vector < typename IterT::value_type::const_iterator>, vector<typename IterT::value_type::const_iterator> >
  get_begins_ends(IterT
  iters_begin,
  IterT iters_end
  ) {
  typedef typename IterT::value_type T;
  typedef vector<typename T::const_iterator> VecT;
  int N = iters_end - iters_begin;
  std::pair<VecT, VecT> result;
  result.first.
  resize(N);
  result.second.
  resize(N);
  for (
  int i = 0;
  i<N;
  i++) {
  result.first[i] = iters_begin[i].

  begin();

  result.second[i] = iters_begin[i].

  end();
}
return
result;
}


template<typename Real, size_t _dim>
class RegularGridInterpolator {

  using GridType = std::vector<std::vector<Real>>;
  using NDArray = boost::multi_array<Real, _dim>;

 public:

  RegularGridInterpolator();

  RegularGridInterpolator(const RegularGridInterpolator &other);

  RegularGridInterpolator &operator=(const RegularGridInterpolator &other);

  void AddCoord(const std::vector<Real> &vec);

  void AddVar(const NDArray &ndarray);

  std::vector<Real> GetCoord(unsigned int idx) const;

  Real GetCoordMinVal(unsigned int idx) const;

  Real GetCoordMaxVal(unsigned int idx) const;

  size_t GetDimension() const;

  Real Interp(const std::array<Real, _dim> &point, bool bound_check=false) const;

  Real EvalAtIndice(const std::array<size_t, _dim> &indice) const;

  bool IsValid() const;

 private:

  void _BuildInterpolator();


 private:
  GridType m_grid_list;
  NDArray m_data;
  std::unique_ptr<InterpMultilinear<_dim, Real, false>> c_interpolator;

  std::array<std::pair<Real, Real>, _dim> m_limits;

};


/*
 * Implementations
 */

template<typename Real, size_t _dim>
RegularGridInterpolator<Real, _dim>::RegularGridInterpolator() {
  m_grid_list.reserve(_dim);
}

template<typename Real, size_t _dim>
RegularGridInterpolator<Real, _dim>::RegularGridInterpolator(const RegularGridInterpolator &other) {
  *this = other;
}

template<typename Real, size_t _dim>
RegularGridInterpolator<Real, _dim> &
RegularGridInterpolator<Real, _dim>::operator=(const RegularGridInterpolator &other) {

  m_grid_list = other.m_grid_list;

  std::vector<size_t> ex;
  const size_t *shape = other.m_data.shape();
  ex.assign(shape, shape + other.m_data.num_dimensions());
  m_data.resize(ex);
  m_data = other.m_data;

  _BuildInterpolator();
  return *this;
}

template<typename Real, size_t _dim>
void RegularGridInterpolator<Real, _dim>::AddCoord(const std::vector<Real> &vec) {
  if (m_grid_list.size() == _dim) {
    std::cerr << "Attempting to add too much dimension to the interpolator" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Check if coord is monoton increasing
  for (size_t i = 1; i < vec.size(); ++i) {
    if (vec[i] <= vec[i-1]) {
      throw std::runtime_error("RegularGridInterpolator only accepts stricly monotonically increasing vectors as grid");
    }
  }

  m_grid_list.push_back(vec);
}

template<typename Real, size_t _dim>
void RegularGridInterpolator<Real, _dim>::AddVar(const NDArray &ndarray) {

  // Verifying that every dimension has been filled
  if (m_grid_list.size() != _dim) {
    std::cerr << "Grid is incompletely defined" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Copying the multi_array
  std::vector<size_t> ex;
  const size_t *shape = ndarray.shape();
  ex.assign(shape, shape + ndarray.num_dimensions());
  m_data.resize(ex);
  m_data = ndarray;

  _BuildInterpolator();

}

template<typename Real, size_t _dim>
void RegularGridInterpolator<Real, _dim>::_BuildInterpolator() {

  // Array of sizes
  std::array<int, _dim> sizes;
  for (int i = 0; i < _dim; i++) {
    sizes[i] = m_grid_list[i].size();
  }

  for (size_t dim=0; dim<_dim; ++dim) {
    auto vec = m_grid_list[dim];
    m_limits[dim].first = m_grid_list[dim][0];
    m_limits[dim].second = m_grid_list[dim][m_grid_list[dim].size()-1];
  }


  auto begins_ends = get_begins_ends(m_grid_list.begin(),
                                     m_grid_list.end()); // TODO: voir a remplacer cet appel copier coller de la doc...
  c_interpolator = std::make_unique<InterpMultilinear<_dim, Real, false>>(begins_ends.first.begin(),
                                                                          sizes,
                                                                          m_data.data(),
                                                                          m_data.data() + m_data.num_elements());
}

template<typename Real, size_t _dim>
Real RegularGridInterpolator<Real, _dim>::Interp(const std::array<Real, _dim> &point, bool bound_check) const {
  // Bound checking
  if (bound_check) {
    for (size_t dim=0; dim<_dim; ++dim) {
      auto coord = point[dim];
      if (coord < m_limits[dim].first || coord > m_limits[dim].second) {
        throw std::runtime_error("RegularGridInterpolator::Interp called with out of bounds value");
      }
    }
  }

  return c_interpolator->interp(point);
}

  template<typename Real, size_t _dim>
  Real RegularGridInterpolator<Real, _dim>::EvalAtIndice(const std::array<size_t, _dim> &indice) const {
    return m_data(indice);
  }

template<typename Real, size_t _dim>
size_t RegularGridInterpolator<Real, _dim>::GetDimension() const {
  return _dim;
}

template<typename Real, size_t _dim>
bool RegularGridInterpolator<Real, _dim>::IsValid() const {
  return c_interpolator != nullptr;
}

template<typename Real, size_t _dim>
std::vector<Real> RegularGridInterpolator<Real, _dim>::GetCoord(unsigned int idx) const {
  if (idx >= _dim) {
    std::cerr << "In RegularGridInterpolator, requesting a coordinate with index that is out of dimension" << std::endl;
    exit(EXIT_FAILURE);
  }
  return m_grid_list[idx];
}

template<typename Real, size_t _dim>
Real RegularGridInterpolator<Real, _dim>::GetCoordMaxVal(unsigned int idx) const {
  if (idx >= _dim) {
    std::cerr << "In RegularGridInterpolator, requesting a coordinate with index that is out of dimension" << std::endl;
    exit(EXIT_FAILURE);
  }
  return m_limits[idx].second;
}

template<typename Real, size_t _dim>
Real RegularGridInterpolator<Real, _dim>::GetCoordMinVal(unsigned int idx) const {
  if (idx >= _dim) {
    std::cerr << "In RegularGridInterpolator, requesting a coordinate with index that is out of dimension" << std::endl;
    exit(EXIT_FAILURE);
  }
  return m_limits[idx].first;
}

}  // end namespace mathutils

#endif //MATHUTILS_MULTILINEARINTERPOLATOR_H
