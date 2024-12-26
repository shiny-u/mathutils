//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_LOOKUPTABLE1D_H
#define MATHUTILS_LOOKUPTABLE1D_H

#include <unordered_map>
#include <iostream>
#include "Interp1d.h"

namespace mathutils {


    // ****************************************************************************************** //
    //                      Class LookupTable1D
    // ****************************************************************************************** //
  template<class Xscalar=double, class Yscalar=double>
  class LookupTable1D {

   protected:
    INTERP_METHOD interp_method;

    std::shared_ptr<std::vector<Xscalar>> Xcoords;
    std::unordered_map<std::string, unsigned long> assoc;
    std::vector<std::shared_ptr<std::vector<Yscalar>>> Ydata;
    std::vector<std::unique_ptr<Interp1d<Xscalar, Yscalar>>>
        interpolators;

   public:
    LookupTable1D(INTERP_METHOD interp_method);

    ~LookupTable1D() {};

    // Le Code qui suit permet de regler le pb de l'utilisation des unique_ptr dans les interpolateurs !!!
    // Voir:
    // https://katyscode.wordpress.com/2012/10/04/c11-using-stdunique_ptr-as-a-class-member-initialization-move-semantics-and-custom-deleters/
    // FIXME: mieux comprendre pourquoi cela fonctionne !!!
    ///Deleting the inplace copy constructor
    LookupTable1D(const LookupTable1D &) = delete;

    /// Deleting the copy operator
    LookupTable1D &operator=(const LookupTable1D &) = delete;

    /// Move constructor allowing to move the interpolators vector storing unique_ptr instances
    LookupTable1D &operator=(LookupTable1D &&table) {
      if (this != &table) {
        Xcoords = std::move(table.Xcoords);
        assoc = std::move(table.assoc);
        Ydata = std::move(table.Ydata);
        interpolators = std::move(table.interpolators);
      }
      return *this;
    }

    /// Move operator allowing to move the interpolators vector storing unique_ptr instances
    LookupTable1D(LookupTable1D &&table)
        : Xcoords(std::move(table.Xcoords)),
          assoc(std::move(table.assoc)),
          Ydata(std::move(table.Ydata)),
          interpolators(std::move(table.interpolators)) {};

    /// Set the X vector of the lookup table
    void SetX(const std::vector<Xscalar> X);

    std::vector<Xscalar> GetX() const { return *Xcoords.get(); }

    Xscalar GetXmin(const std::string &name) const { return interpolators.at(GetIndex(name))->GetXmin(); }

    Xscalar GetXmax(const std::string &name) const { return interpolators.at(GetIndex(name))->GetXmax(); }

    /// Get the number of series
    unsigned long GetNbSeries() const { return Ydata.size(); }

    /// Get the number of samples
    unsigned long GetNbSample() const { return Xcoords->size(); }

    /// Add a Serie to the LUT
    bool AddY(const std::string name, const std::vector<Yscalar> Y);

    /// Evaluates the LUT giving the key of the serie and a value
    Yscalar Eval(const std::string name, const Xscalar x) const;

    /// Evaluates the Grad giving the key of the serie and a value
    Yscalar GradEval(const std::string name, const Xscalar x) const;

    /// Evaluates the LUT giving the key of the serie and a vector of values
    std::vector<Yscalar> Eval(const std::string name, const std::vector<Xscalar> &xvect) const;

    /// Evaluates the LUT giving the key of the serie and a vector of values
    std::vector<Yscalar> GradEval(const std::string name, const std::vector<Xscalar> &xvect) const;

    /// Returns true if the serie exists in the table
    bool Exists(const std::string &name) const;

    /// Evaluates the LUT giving a value
    template<class Tx, class Ty>
    std::unordered_map<std::string, Ty> Eval(const Tx x) const;

   private:
    /// Get the index of the series from its name
    inline unsigned long GetIndex(const std::string name) const;
  };


  template<class Xscalar, class Yscalar>
  LookupTable1D<Xscalar, Yscalar>::LookupTable1D(INTERP_METHOD interp_method):
      interp_method(interp_method) {

  }

  template<class Xscalar, class Yscalar>
  void LookupTable1D<Xscalar, Yscalar>::SetX(const std::vector<Xscalar> X) {
    Xcoords = std::make_shared<std::vector<Xscalar>>(X);
  }

  template<class Xscalar, class Yscalar>
  bool LookupTable1D<Xscalar, Yscalar>::AddY(const std::string name, const std::vector<Yscalar> Y) {
    // TODO: verifier qu'on a le meme nombre d'elt que dans X...
    // Get the future position of the new Serie
    auto i = GetNbSeries();

    // Trying to insert the new name/index into the association map
    auto res_pair = assoc.insert(std::pair<std::string, unsigned long>(name, i));

    if (!res_pair.second) { // Name already exists into the map
      std::cout << "Data have not been added to the LUT" << std::endl;

    } else {
      // We can add data
      auto Y_shared = std::make_shared<std::vector<Yscalar>>(Y);
      Ydata.push_back(Y_shared);

      // Building the interpolator based on the global interpolation method of the LUT
      auto interp_ptr = Interp1d<Xscalar, Yscalar>::MakeInterp1d(interp_method);
      // FIXME: ICI, on ne recupere pas comme voulu un pointeur vers un objet FrInterp1dLinear
      // Du coup, la methode Initialize appelee apres n'est que celle de

      // Initializing the interpolator
      interp_ptr->Initialize(Xcoords, Y_shared);

      auto interp_unique = std::unique_ptr<Interp1d<Xscalar, Yscalar>>(interp_ptr);

      interpolators.push_back(std::move(interp_unique));

    }

    return res_pair.second;
  }

  template<class Xscalar, class Yscalar>
  Yscalar LookupTable1D<Xscalar, Yscalar>::Eval(const std::string name, const Xscalar x) const {
    return (*interpolators.at(GetIndex(name)))(x);
  }

  template<class Xscalar, class Yscalar>
  Yscalar LookupTable1D<Xscalar, Yscalar>::GradEval(const std::string name, const Xscalar x) const {
    return (*interpolators.at(GetIndex(name))).GradEval(x);
  }

  template<class Xscalar, class Yscalar>
  std::vector<Yscalar>
  LookupTable1D<Xscalar, Yscalar>::Eval(const std::string name, const std::vector<Xscalar> &xvect) const {
    return interpolators.at(GetIndex(name))->Eval(xvect);
  }

  template<class Xscalar, class Yscalar>
  std::vector<Yscalar>
  LookupTable1D<Xscalar, Yscalar>::GradEval(const std::string name, const std::vector<Xscalar> &xvect) const {
    return interpolators.at(GetIndex(name))->GradEval(xvect);
  }

  template<class Xscalar, class Yscalar>
  bool LookupTable1D<Xscalar, Yscalar>::Exists(const std::string &name) const {
    return assoc.find(name) != assoc.end();
  }

  template<class Xscalar, class Yscalar>
  template<class Tx, class Ty>
  std::unordered_map<std::string, Ty> LookupTable1D<Xscalar, Yscalar>::Eval(const Tx x) const {

    std::unordered_map<std::string, Ty> out;
    out.reserve(GetNbSeries());

    std::string name;
    unsigned long idx;

    for (auto elt: assoc) {
      name = elt.first;
      idx = elt.second;

      out[name] = interpolators[idx]->Eval(x);

    }
    return out;

  }

  template<class Xscalar, class Yscalar>
  unsigned long LookupTable1D<Xscalar, Yscalar>::GetIndex(const std::string name) const {
    return assoc.at(name);
  }



    // ****************************************************************************************** //
    //                      Class LookupTable1D_DataContained
    // ****************************************************************************************** //
    template<class Xscalar=double, class Yscalar=double>
    class LookupTable1D_DataContained {
    protected:
        INTERP_METHOD interp_method;
        std::vector<Xscalar> Xcoords;
        std::unordered_map<std::string, unsigned long> assoc;
        std::vector<std::vector<Yscalar>> Ydata;
        std::vector<std::shared_ptr<Interp1d<Xscalar, Yscalar>>> interpolators;

    public:
        LookupTable1D_DataContained(INTERP_METHOD interp_method) : interp_method(interp_method) {}

        // Copy constructor
        LookupTable1D_DataContained(const LookupTable1D_DataContained& other) :
                interp_method(other.interp_method),
                Xcoords(other.Xcoords),
                assoc(other.assoc),
                Ydata(other.Ydata)
        {
            // Recreate interpolators based on the interp_method
            int iY = 0;
            for (const auto& interp : other.interpolators)
            {
                std::shared_ptr<Interp1d<Xscalar, Yscalar>> new_interp;
                switch (interp_method) {
                    case LINEAR:
                        new_interp = std::make_shared<Interp1dLinear<Xscalar, Yscalar>>();
                        break;
                    case LINEAR_EXTRAPOLATE:
                        new_interp = std::make_shared<Interp1dLinearExtrapolate<Xscalar, Yscalar>>();
                        break;
                    case LINEAR_SATURATE:
                        new_interp = std::make_shared<Interp1dLinearSaturate<Xscalar, Yscalar>>();
                        break;
                    default:
                        // Handle default case
                        break;
                }
                if (new_interp)
                {
                    auto Y_shared = std::make_shared<std::vector<Yscalar>>(Ydata[iY]);
                    auto Xcoords_shared = std::make_shared<std::vector<Xscalar>>(Xcoords);
                    new_interp->Initialize(Xcoords_shared, Y_shared);
                    interpolators.push_back(new_interp);
                }
                iY++;
            }
        }

        // Copy assignment operator
        LookupTable1D_DataContained& operator=(const LookupTable1D_DataContained& other) {
            if (this != &other) {
                interp_method = other.interp_method;
                Xcoords = other.Xcoords;
                assoc = other.assoc;
                Ydata = other.Ydata;
                interpolators.clear(); // Clear existing interpolators

                // Recreate interpolators based on the interp_method
                int iY = 0;
                for (const auto& interp : other.interpolators)
                {
                    std::shared_ptr<Interp1d<Xscalar, Yscalar>> new_interp;
                    switch (interp_method) {
                        case LINEAR:
                            new_interp = std::make_shared<Interp1dLinear<Xscalar, Yscalar>>();
                            break;
                        case LINEAR_EXTRAPOLATE:
                            new_interp = std::make_shared<Interp1dLinearExtrapolate<Xscalar, Yscalar>>();
                            break;
                        case LINEAR_SATURATE:
                            new_interp = std::make_shared<Interp1dLinearSaturate<Xscalar, Yscalar>>();
                            break;
                        default:
                            // Handle default case
                            break;
                    }
                    if (new_interp)
                    {
                        auto Y_shared = std::make_shared<std::vector<Yscalar>>(Ydata[iY]);
                        auto Xcoords_shared = std::make_shared<std::vector<Xscalar>>(Xcoords);
                        new_interp->Initialize(Xcoords_shared, Y_shared);
                        interpolators.push_back(new_interp);
                    }
                    iY++;
                }
            }
            return *this;
        }

        void SetX(const std::vector<Xscalar> &X) {
            Xcoords = X;
        }

        std::vector<Xscalar> GetX() const {
            return Xcoords;
        }

        INTERP_METHOD get_interp_method() const{
            return interp_method;
        }

        Xscalar GetXmin(const std::string &name) const {
            return interpolators.at(GetIndex(name)).GetXmin();
        }

        Xscalar GetXmax(const std::string &name) const {
            return interpolators.at(GetIndex(name)).GetXmax();
        }

        unsigned long GetNbSeries() const {
            return Ydata.size();
        }

        unsigned long GetNbSample() const {
            return Xcoords.size();
        }

        bool AddY(const std::string &name, const std::vector<Yscalar> &Y) {
            auto i = GetNbSeries();
            auto res_pair = assoc.insert(std::pair<std::string, unsigned long>(name, i));
            if (!res_pair.second) {
                std::cout << "Data have not been added to the LUT" << std::endl;
                return false;
            }
            else
            {
                auto Y_shared = std::make_shared<std::vector<Yscalar>>(Y);
                auto Xcoords_shared = std::make_shared<std::vector<Xscalar>>(Xcoords);
                Ydata.push_back(Y);
                switch (interp_method)
                {
                    case LINEAR:{
                        auto interp = std::make_shared<Interp1dLinear<Xscalar, Yscalar>>();
                        interp->Initialize(Xcoords_shared, Y_shared);
                        interpolators.push_back(interp);
                        break;
                    }
                    case LINEAR_EXTRAPOLATE:
                    {
                        auto interp = std::make_shared<Interp1dLinearExtrapolate<Xscalar, Yscalar>>();
                        interp->Initialize(Xcoords_shared, Y_shared);
                        interpolators.push_back(interp);
                        break;
                    }
                    case LINEAR_SATURATE:
                    {
                        auto interp = std::make_shared<Interp1dLinearSaturate<Xscalar, Yscalar>>();
                        interp->Initialize(Xcoords_shared, Y_shared);
                        interpolators.push_back(interp);
                        break;
                    }
                    default:
                    {
                        // DO NOTHING
                        break;
                    }
                }
                return true;
            }
        }

        Yscalar Eval(const std::string &name, const Xscalar x) const {
            return interpolators.at(GetIndex(name))->Eval(x);
        }

        Yscalar GradEval(const std::string &name, const Xscalar x) const {
            return interpolators.at(GetIndex(name))->GradEval(x);
        }

        std::vector<Yscalar> Eval(const std::string &name, const std::vector<Xscalar> &xvect) const {
            return interpolators.at(GetIndex(name))->val(xvect);
        }

        std::vector<Yscalar> GradEval(const std::string &name, const std::vector<Xscalar> &xvect) const {
            return interpolators.at(GetIndex(name))->GradEval(xvect);
        }

        bool Exists(const std::string &name) const {
            return assoc.find(name) != assoc.end();
        }

        template<class Tx, class Ty>
        std::unordered_map<std::string, Ty> Eval(const Tx x) const {
            std::unordered_map<std::string, Ty> out;
            out.reserve(GetNbSeries());
            for (const auto &elt: assoc) {
                const std::string &name = elt.first;
                unsigned long idx = elt.second;
                out[name] = interpolators[idx](x);
            }
            return out;
        }

    private:
        unsigned long GetIndex(const std::string &name) const {
            return assoc.at(name);
        }
    };
}  // end namespace mathutils



#endif //MATHUTILS_LOOKUPTABLE1D_H
