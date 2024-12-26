//
// Created by frongere on 29/11/17.
//

#ifndef MATHUTILS_INTERP2D_H
#define MATHUTILS_INTERP2D_H

#include "Interp.h"

// TODO: rendre capable la methode d'inteprolation 2D d'interpoler les nombres complexes !!

namespace mathutils {

    template <class Real>
    class Interp2d {

    protected:
        std::shared_ptr<const std::vector<Real>> xcoord;
        std::shared_ptr<const std::vector<Real>> ycoord;
        std::shared_ptr<const std::vector<Real>> mdata;
        unsigned long nx = 0;
        Real xmin;
        Real xmax;
        unsigned long ny = 0;
        Real ymin;
        Real ymax;
        unsigned long ndata = 0;

    public:
//        Interp2d() {};
//        ~Interp2d() {};

        virtual void Initialize(std::shared_ptr<const std::vector<Real>> x,
                                std::shared_ptr<const std::vector<Real>> y,
                                std::shared_ptr<const std::vector<Real>> data);

        virtual Real Eval(const Real x, const Real y) const = 0;

        //virtual Real Eval(const std::vector<Real> coord) const = 0;

        virtual std::vector<Real> Eval(const std::vector<std::vector<Real>> vcoord) const = 0;

        Real operator() (const Real x, const Real y) const { return Eval(x, y); }

        std::vector<Real> operator() (const std::vector<std::vector<Real>> vcoord) const {
            return Eval(vcoord);
        }

        static Interp2d<Real>* MakeInterp2d(INTERP_METHOD method);

        Real GetXmin() const { return xmin; }

        Real GetXmax() const { return xmax; }

        Real GetYmin() const { return ymin; }

        Real GetYmax() const { return ymax; }

    };


    template <class Real>
    void Interp2d<Real>::Initialize(std::shared_ptr<const std::vector<Real>> x,
                                    std::shared_ptr<const std::vector<Real>> y,
                                    std::shared_ptr<const std::vector<Real>> data) {

        //assert( x->size() * y->size() == data->size());
        //assert( y->size() == data->at(0).size());
        assert( std::is_sorted(x->begin(), x->end()) );
        assert( std::is_sorted(y->begin(), y->end()) );

        nx = x->size();
        xmin = x->at(0);
        xmax = x->at(nx-1);

        ny = y->size();
        ymin = y->at(0);
        ymax = y->at(ny-1);

        ndata = nx*ny;

        xcoord = x;
        ycoord = y;
        mdata = data;

    }

    template <class Real>
    class Interp2dLinear : public Interp2d<Real> {

    private:
        std::vector<Real> a1;
        std::vector<Real> a2;
        std::vector<Real> a3;
        std::vector<Real> a4;

    public:
        void Initialize(const std::shared_ptr<const std::vector<Real>> x,
                        const std::shared_ptr<const std::vector<Real>> y,
                        const std::shared_ptr<const std::vector<Real>> data) override;

        Real Eval(const Real x, const Real y) const;

        //Real Eval(const std::vector<Real> coord) const;

        std::vector<Real> Eval(const std::vector<std::vector<Real>> vcoord) const;

    };

    template <class Real>
    void Interp2dLinear<Real>::Initialize(const std::shared_ptr<const std::vector<Real>> x,
                                          const std::shared_ptr<const std::vector<Real>> y,
                                          const std::shared_ptr<const std::vector<Real>> data) {

        Interp2d<Real>::Initialize(x, y, data);

        a1.reserve(this->ndata);
        a2.reserve(this->ndata);
        a3.reserve(this->ndata);
        a4.reserve(this->ndata);

        auto ny = this->ny;

        Real x1, x2, y1, y2, inv_dx, inv_dy;
        Real data_11, data_12, data_21, data_22;
        std::shared_ptr<std::vector<Real>> data_1, data_2;

        for (unsigned int i = 1; i < this->nx; i++) {

            x1 = this->xcoord->at(i-1);
            x2 = this->xcoord->at(i);

            inv_dx = 1. / (x2 - x1);

            //data_1 = std::make_shared<std::vector<Real>>(data->at(i-1));
            //data_2 = std::make_shared<std::vector<Real>>(data->at(i));

            for (unsigned int j = 1; j < this->ny; j++) {

                y1 = this->ycoord->at(j-1);
                y2 = this->ycoord->at(j);

                inv_dy = 1./ (y2 - y1);

                data_11 = data->at(ny*(i-1)+(j-1));
                data_12 = data->at(ny*(i-1)+j);
                data_21 = data->at(ny*i+(j-1));
                data_22 = data->at(ny*i+j);

                a1.push_back( (data_11*x2*y2 - data_12*x2*y1 - data_21*x1*y2 + data_22*x1*y1) *
                              inv_dx * inv_dy );

                a2.push_back( (-data_11*y2 + data_12*y1 + data_21*y2 - data_22*y1) *
                              inv_dx * inv_dy );

                a3.push_back( (-data_11*x2 + data_12*x2 + data_21*x1 - data_22*x1) *
                              inv_dx * inv_dy );

                a4.push_back( (data_11 - data_12 - data_21 + data_22) * inv_dx * inv_dy);

            }
        }

    }

    template <class Real>
    Real Interp2dLinear<Real>::Eval(const Real x, const Real y) const {

        assert (x >= this->xmin && x <= this->xmax);
        assert (y >= this->ymin && y <= this->ymax);

        auto upper_x = std::lower_bound(this->xcoord->begin(), this->xcoord->end(), x);
        auto ix = std::distance(this->xcoord->begin(), upper_x);

        auto upper_y = std::lower_bound(this->ycoord->begin(), this->ycoord->end(), y);
        auto iy = std::distance(this->ycoord->begin(), upper_y);

        if (ix == 0) ix = 1;
        if (iy == 0) iy = 1;

        auto idata = iy + (ix-1)*(this->ny-1) - 1;

        Real a1_ = a1.at(idata);
        Real a2_ = a2.at(idata);
        Real a3_ = a3.at(idata);
        Real a4_ = a4.at(idata);

        return a1_ + a2_*x + a3_*y + a4_*x*y;
    }

    //template <class Real>
    //Real Interp2dLinear<Real>::Eval(const std::vector<Real> coord) const {
    //    return Eval(coord[0], coord[1]);
    //}

    template <class Real>
    std::vector<Real> Interp2dLinear<Real>::Eval(const std::vector<std::vector<Real>> vcoord ) const {

        auto n = vcoord.size();

        std::vector<Real> out;
        out.reserve(n);

        for (unsigned long i = 0; i < n; i++) {
            out.push_back(Eval(vcoord[i][0], vcoord[i][1]));
        }

        return out;
    }


    template <class Real>
    Interp2d<Real>* Interp2d<Real>::MakeInterp2d(INTERP_METHOD method) {
        switch (method) {
            case LINEAR:
                return new Interp2dLinear<Real>;
            default:
                throw("2D INTERPOLATION METHOD DOES NOT EXIST");
        }
    }

}  // end namespace mathutils

#endif //MATHUTILS_INTERP2D_H
