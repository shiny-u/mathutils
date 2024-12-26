//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_INTERP1D_H
#define MATHUTILS_INTERP1D_H

#include "Interp.h"
#include "Spline.h"
#include <iostream>

namespace mathutils {


  template<class XReal, class YScalar>
  class Interp1d {

   protected:
    std::shared_ptr<const std::vector<XReal>> xcoord;
    std::shared_ptr<const std::vector<YScalar>> yval;
    unsigned long ndata = 0;
    XReal xmin;
    XReal xmax;

   protected:

    virtual YScalar PrivEval(XReal x) const = 0;

    virtual YScalar PrivGradEval(XReal x) const = 0;

   public:

    virtual ~Interp1d() = default;

    YScalar Eval(XReal x) const {return PrivEval(x);};

    YScalar GradEval(XReal x) const {return PrivGradEval(x);};

    std::vector<YScalar> Eval(const std::vector<XReal> &xvector) const {

      auto n = xvector.size();

      std::vector<YScalar> out;
      out.reserve(n);

      for (int i = 0; i < n; i++) {
        out.push_back(Eval(xvector[i]));
      }
      return out;
    }


    std::vector<YScalar> GradEval(const std::vector<XReal> &xvector) const {
      auto n = xvector.size();
      std::vector<YScalar> out;
      out.reserve(n);
      for (int i = 0; i < n; i++) {
         out.push_back(GradEval(xvector[i]));
      }
      return out;
    }

   public:
    // TODO: voir a separer l'implementation et la mettre en fin de fichier (pas directement dans le corps de la classe)

    virtual void Initialize(std::shared_ptr<const std::vector<XReal>> x,
                            std::shared_ptr<const std::vector<YScalar>> y);

    YScalar operator()(const XReal x) const {
      return Eval(x);
    }

    std::vector<YScalar> operator()(const std::vector<XReal> &xvector) const {
      return Eval(xvector);
    }

    static Interp1d<XReal, YScalar> *MakeInterp1d(INTERP_METHOD method);

    XReal GetXmin() const { return xmin; }

    XReal GetXmax() const { return xmax; }

  };

  template<class XReal, class YScalar>
  void Interp1d<XReal, YScalar>::Initialize(std::shared_ptr<const std::vector<XReal>> x,
                                            std::shared_ptr<const std::vector<YScalar>> y) {

    assert(x->size() == y->size());
    assert (std::is_sorted(x->begin(), x->end()));

    ndata = x->size();
    xmin = x->at(0);
    xmax = x->at(ndata - 1);

    xcoord = x;
    yval = y;

  }


  template<class XReal, class YScalar>
  class Interp1dLinearBase : public Interp1d<XReal, YScalar> {

   private:
    std::vector<YScalar> a;
    std::vector<YScalar> b;

   public:

    void Initialize(std::shared_ptr<const std::vector<XReal>> x,
                    std::shared_ptr<const std::vector<YScalar>> y) override;

   protected:
    YScalar PrivEval(XReal x) const override;

    YScalar PrivGradEval(XReal x) const override;
//    std::vector<YScalar> Eval(const std::vector<XReal> &xvector) const override;

  };


  template<class XReal, class YScalar>
  class Interp1dLinear : public Interp1dLinearBase<XReal, YScalar> {

   protected:

    YScalar PrivEval(XReal x) const override;
    YScalar PrivGradEval(XReal x) const override;

  };


  template<class XReal, class YScalar>
  class Interp1dLinearExtrapolate : public Interp1dLinearBase<XReal, YScalar> {

   protected:

    YScalar PrivEval(XReal x) const override;

  };


  template<class XReal, class YScalar>
  class Interp1dLinearSaturate : public Interp1dLinearBase<XReal, YScalar> {

   protected:

    YScalar PrivEval(XReal x) const override;
    YScalar PrivGradEval(XReal x) const override;
  };


  template<class XReal, class YScalar>
  class Interp1dBSpline : public Interp1d<XReal, YScalar> {
   private:
//        Spline m_spline;

   public:
    // TODO: terminer
  };


  template<class XReal, class YScalar>
  void Interp1dLinearBase<XReal, YScalar>::Initialize(const std::shared_ptr<const std::vector<XReal>> x,
                                                      const std::shared_ptr<const std::vector<YScalar>> y) {

    Interp1d<XReal, YScalar>::Initialize(x, y);

    a.reserve(this->ndata);
    b.reserve(this->ndata);

    XReal xi, xii, xii_m_xi;
    YScalar yi, yii;
    for (unsigned int i = 1; i < this->ndata; ++i) {

      xi = this->xcoord->at(i - 1);
      xii = this->xcoord->at(i);

      yi = this->yval->at(i - 1);
      yii = this->yval->at(i);

      xii_m_xi = xii - xi;

      a.push_back((yii - yi) / xii_m_xi);
      b.push_back((yi * xii - xi * yii) / xii_m_xi);

    }
  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinearBase<XReal, YScalar>::PrivEval(const XReal x) const {
    // TODO: il faut que le type de retour soit compatible avec real et complex !!!

    // First, binary search on the x coords
    auto upper = std::lower_bound(this->xcoord->begin(), this->xcoord->end(), x);
    auto index = std::distance(this->xcoord->begin(), upper);

    if (index == 0) index = 1;  // Bug fix for x == xmin
    if (index > a.size()) index = a.size();

    YScalar a_ = a.at(index - 1);
    YScalar b_ = b.at(index - 1);

    return a_ * x + b_;
  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinearBase<XReal, YScalar>::PrivGradEval(const XReal x) const {

    // First, binary search on the x coords
    auto upper = std::lower_bound(this->xcoord->begin(), this->xcoord->end(), x);
    auto index = std::distance(this->xcoord->begin(), upper);

    // Bug fix for x == xmin
    if (index == 0)
        index = 1;

    if (index > a.size())
        index = a.size();

    return a.at(index - 1);
  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinear<XReal, YScalar>::PrivEval(const XReal x) const {

    if (x < this->xmin || x > this->xmax) {
      std::cerr << "Interpolation evaluated for value " << x << ", outside of the range : [" << this->xmin << ", "
                << this->xmax << "]" << std::endl;
      exit(1);
    }

    return Interp1dLinearBase<XReal, YScalar>::PrivEval(x);

  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinear<XReal, YScalar>::PrivGradEval(const XReal x) const {
    if (x < this->xmin || x > this->xmax) {
        std::cerr << "Interpolation evaluated for value " << x << ", outside of the range : [" << this->xmin << ", "
                  << this->xmax << "]" << std::endl;
        exit(1);
    }

    return Interp1dLinearBase<XReal, YScalar>::PrivGradEval(x);
  }
  template<class XReal, class YScalar>
  YScalar Interp1dLinearExtrapolate<XReal, YScalar>::PrivEval(const XReal x) const {
    return Interp1dLinearBase<XReal, YScalar>::PrivEval(x);
  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinearSaturate<XReal, YScalar>::PrivEval(const XReal x) const {

    auto x_tmp = x;
    if (x < this->xmin || x > this->xmax) {
      x < this->xmin ? x_tmp = this->xmin : x_tmp = this->xmax;
    }
    return Interp1dLinearBase<XReal, YScalar>::PrivEval(x_tmp);

  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinearSaturate<XReal, YScalar>::PrivGradEval(const XReal x) const {
      if (x < this->xmin)
          return YScalar{};
      else if (x > this->xmax)
          return YScalar{};
      else
          return Interp1dLinearBase<XReal, YScalar>::PrivGradEval(x);
   }
//
//  template<class XReal, class YScalar>
//  std::vector<YScalar> Interp1d<XReal, YScalar>::Eval(const std::vector<XReal> &xvector) const

  /// Factory method to create 1D interpolation classes
  template<class XReal, class YScalar>
  Interp1d<XReal, YScalar> *Interp1d<XReal, YScalar>::MakeInterp1d(INTERP_METHOD method) {
    Interp1d<XReal, YScalar> *interp1D;
    switch (method) {
      case LINEAR:
        interp1D = new Interp1dLinear<XReal, YScalar>;
        break;
      case LINEAR_EXTRAPOLATE:
        interp1D = new Interp1dLinearExtrapolate<XReal, YScalar>;
        break;
      case LINEAR_SATURATE:
        interp1D = new Interp1dLinearSaturate<XReal, YScalar>;
        break;
      case BSPLINE:
//        return new Interp1dBSpline<XReal, YScalar>;
        std::cerr << "no B-Slpine interpolator available yet" << std::endl;
        exit(1);
      default:
        std::cerr << "1D INTERPOLATION METHOD DOES NOT EXIST" << std::endl;
        exit(1);
    }
    return interp1D;
  }

}  // end namespace mathutils



#endif //MATHUTILS_INTERP1D_H
