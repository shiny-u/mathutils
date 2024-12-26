//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_LOOKUPTABLE2D_H
#define MATHUTILS_LOOKUPTABLE2D_H

#include <vector>

#include <unordered_map>

#include <iostream>

#include "Interp2d.h"

namespace mathutils {

//    template <class Real>
//    class Interp2d<Real>;

    template <class Real=double>
    class LookupTable2d {

    protected:
        INTERP_METHOD interp_method = LINEAR;

        std::shared_ptr<std::vector<Real>> Xcoords;
        std::shared_ptr<std::vector<Real>> Ycoords;
        std::vector<std::shared_ptr<std::vector<Real>>> MData;
        std::unordered_map<std::string, unsigned long> assoc;
        std::vector<std::unique_ptr<Interp2d<Real>>> interpolators;

    public:

        LookupTable2d() {};
        ~LookupTable2d() {};

        /// Deleting the inplace copy constructor
        //LookupTable2d(const LookupTable2d &) = delete;

        /// Deleting the copy operator
        //LookupTable2d & operator = (const LookupTable2d &) = delete;

        /// Set interpolation method
        void SetInterpolationMethod(INTERP_METHOD method);

        /// Set the X vector of the lookup table
        void SetX(const std::vector<Real> X);

        std::vector<Real> GetX() const { return *Xcoords.get(); }

        /// Set the Y vector of the lookup table
        void SetY(const std::vector<Real> Y);

        std::vector<Real> GetY() const { return *Ycoords.get(); }

        /// Get the number of series
        unsigned long GetNbSeries() const { return MData.size(); }

        /// Get the number of samples
        unsigned long GetNbSample() const;

        /// Add a serie to the LUT
        bool AddData(std::string name, std::vector<Real> Data);

        /// Evaluates the LUT giving the key of the serie and the coordinate of a point
        Real Eval(const std::string name, const Real x, const Real y) const;
        Real Eval(const std::string name, const std::vector<Real> coord) const;

        /// Evaluates the LUT giving the key of the series and a vector of coordinate
        std::vector<Real> Eval(const std::string name, const std::vector<std::vector<Real>>& vcoord) const;

        bool HasSerie(const std::string name) const;

    private:

        /// Get the index fo the series from its name
        inline unsigned long GetIndex(const std::string name) const;

    };

    template <class Real>
    void LookupTable2d<Real>::SetInterpolationMethod(INTERP_METHOD method) {

        if (GetNbSeries() > 0) {
            std::cout << "TODO : reinitialiser les interpolateurs" << std::endl;
        }
        interp_method = method;
    }

    template <class Real>
    unsigned long LookupTable2d<Real>::GetNbSample() const {
        return Xcoords->size() * Ycoords->size();
    }

    template <class Real>
    void LookupTable2d<Real>::SetX(const std::vector<Real> X) {
        Xcoords = std::make_shared<std::vector<Real>>(X);
    }

    template <class Real>
    void LookupTable2d<Real>::SetY(const std::vector<Real> Y) {
        Ycoords = std::make_shared<std::vector<Real>>(Y);
    }

    template <class Real>
    bool LookupTable2d<Real>::AddData(const std::string name, const std::vector<Real> Data) {

        auto i = GetNbSeries();

        auto res_pair = assoc.insert(std::pair<std::string, unsigned long>(name, i));

        if (!res_pair.second) {
            std::cout << "Data have not been added to the LUT" << std::endl;
        } else {

            auto Data_shared = std::make_shared<std::vector<Real>>(Data);
            MData.push_back(Data_shared);

            auto interp_ptr = Interp2d<Real>::MakeInterp2d(interp_method);

            interp_ptr->Initialize(Xcoords, Ycoords, Data_shared);

            auto interp_unique = std::unique_ptr<Interp2d< Real>>(interp_ptr);

            interpolators.push_back(std::move(interp_unique));
        }

        return res_pair.second;
    }

    template <class Real>
    Real LookupTable2d<Real>::Eval(const std::string name, const Real x, const Real y) const {
        return interpolators.at(GetIndex(name))->Eval(x, y);
    }

    template <class Real>
    Real LookupTable2d<Real>::Eval(const std::string name, const std::vector<Real> coord) const {
        return interpolators.at(GetIndex(name))->Eval(coord);
    }

    template <class Real>
    std::vector<Real> LookupTable2d<Real>::Eval(const std::string name, const std::vector<std::vector<Real>>& vcoord) const {
        return interpolators.at(GetIndex(name))->Eval(vcoord);
    }

    template <class Real>
    bool LookupTable2d<Real>::HasSerie(const std::string name) const {

        if (assoc.find(name) == assoc.end()) {
            return false;
        } else {
            return true;
        }
    }

    template <class Real>
    unsigned long LookupTable2d<Real>::GetIndex(const std::string name) const {
        return assoc.at(name);
    }

}  // end namespace mathutils

#endif //MATHUTILS_LOOKUPTABLE2D_H
