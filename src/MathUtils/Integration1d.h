//
// Created by frongere on 15/12/17.
//

#ifndef MATHUTILS_INTEGRATION1D_H
#define MATHUTILS_INTEGRATION1D_H

#include <vector>

namespace mathutils {

    enum INTEGRATION_METHOD {
        TRAPEZOIDAL,
        NEWTON_COTES
    };


    template <class Scalar>
    class Integrate1d {

        enum INTEGRATION_DATA {
            VECTOR,
            FUNCTION
        };

    protected:
        Scalar m_Xmin = 0.;
        Scalar m_Xmax = 1.;
        unsigned int m_NbPoints;

        std::vector<Scalar> m_Y;
        Scalar (*m_integrand)(Scalar x) = nullptr; // Function to be integrated

        INTEGRATION_METHOD m_IntegrationMethod = TRAPEZOIDAL;
        INTEGRATION_DATA c_DataSource = VECTOR;

        bool c_Computed = false;

        Scalar c_Result = 0.;


    public:
        Integrate1d(Scalar (*F)(Scalar x), Scalar xmin, Scalar xmax, unsigned int nbPoints)
                : m_Xmin(xmin),
                  m_Xmax(xmax),
                  m_NbPoints(nbPoints),
                  c_Computed(false) {

            m_integrand = F;
            c_DataSource = FUNCTION;
        }

        Integrate1d(const std::vector<Scalar>& Y, Scalar xmin, Scalar xmax, unsigned int nbPoints) :
                m_Xmin(xmin),
                m_Xmax(xmax),
                m_NbPoints(nbPoints),
                c_Computed(false) {  // FIXME: pas besoin de mettre nbPoints...

            m_Y = Y;
            c_DataSource = VECTOR;
        }

        void SetIntegrationMethod(INTEGRATION_METHOD method) { m_IntegrationMethod = method; }

        INTEGRATION_METHOD GetIntegrationMethod() const { return m_IntegrationMethod; }

        void SetXmin(Scalar xmin) {
            m_Xmin = xmin;
            c_Computed = false;
        }

        Scalar GetXmin() const { return m_Xmin; }

        void SetXmax(Scalar xmax) {
            m_Xmax = xmax;
            c_Computed = false;
        }

        Scalar GetXmax() const {
            return m_Xmax;
        }

        Scalar Get() {
            if (!c_Computed) {
                Compute();
            }
            return c_Result;
        }

    private:
        void Compute() {

            switch (m_IntegrationMethod) {
                case TRAPEZOIDAL:
                    ComputeTrapz();
                    break;
                case NEWTON_COTES:
                    ComputeNewtonCotes();
                    break;
            }

            c_Computed = true;
        }

        void ComputeTrapz();
        void ComputeNewtonCotes();
    };


    // Trapezoidal rule
    template <class Scalar>
    void Integrate1d<Scalar>::ComputeTrapz() {

        Scalar dx = (m_Xmax-m_Xmin) / Scalar(m_NbPoints-1);

        if (c_DataSource == FUNCTION) {
            m_Y.clear();
            m_Y.swap(m_Y);

            Scalar x(m_Xmin);

            for (unsigned int i=0; i<m_NbPoints; i++) {
                m_Y.push_back((*m_integrand)(x));
                x += dx;
            }
        }

        Scalar sum = 0.0;
        for (unsigned int i=1; i<m_NbPoints-1; i++) {
            sum += m_Y[i];
        }

        c_Result = m_Y[0] + m_Y[m_Y.size()-1] + 2.0*sum;  // TODO: drole de costruction: faire sum + 0.5*(f[0] + f[1], plus simple
        c_Result *= 0.5 * dx;

    }

    template <class Scalar>
    void Integrate1d<Scalar>::ComputeNewtonCotes() {
        // TODO
        throw std::runtime_error("Newton-Cotes 1D integration method is not implemented yet");
    }


    // =================================================================================================================
    // =================================================================================================================
    // INTEGRATION FUNCTIONS
    // =================================================================================================================
    // =================================================================================================================

    template <class Scalar>
    Scalar Trapz(const std::vector<Scalar>& y, Scalar& dx=1.) { // TODO: voir a utiliser cette fonction dans la classe ci-dessus
        // This is the particular case where we have a constant step size between samples
        // Simplification is worth a new implementation

        unsigned long N = y.size();

        Scalar sum = 0.;
        for (unsigned long i = 1; i<N-1; i++) {
            sum += y[i];
        }

        return dx * (sum + 0.5 * (y[0] + y[N-1]));
    }

    template <class Scalar>
    Scalar Trapz(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
        // This implementation is suitable for non-uniform sampling

        unsigned long N = y.size();

        assert(N > 1);
        assert(x.size() == N);

        Scalar dx1 = x[1] - x[0];
        Scalar dxN_1 = x[N-1] - x[N-2];

        Scalar sum = 0.;
        Scalar dxi, dxii;

        dxii = dx1;

        for (unsigned long i=1; i<N-1; i++) {
            dxi = dxii;
            dxii = x[i+1] - x[i];
            sum += y[i] * (dxi + dxii);
        }

        return 0.5 * (sum + y[0]*dx1 + y[N-1]*dxN_1);

    }







}  // end namespace mathutils


#endif //MATHUTILS_INTEGRATION1D_H
