#include <complex>
#include <vector>
#include <iostream>


#define OMEGA_FUN(name) __ ## opr_ww_i1 ## _MOD_ ## name


extern "C" {
    void OMEGA_FUN(new_event)(const double*);
    double OMEGA_FUN(color_sum)(int& flv, int& hel);
    void OMEGA_FUN(init)(const double*, int& scheme);
}

class OmegaWrapper {
    public:
    OmegaWrapper(const std::vector<double>& parameters) : m_parameters(parameters) {}

    double operator()(std::vector<double> lvec, int flv)
    {
        return sqme(lvec.data(), flv);
    }

    double sqme(double* lvec, int flv)
    {
        init_parameters();
        new_event(lvec);

        double res = 0;
        for (int i = 1; i <= 64; i++) {
            res += OMEGA_FUN(color_sum)(flv, i);
        }
        //average over input helicity combinations
        return res / 4.0;
    }


    private:
    const std::vector<double> m_parameters;
    void init_parameters() {
        int scheme = 0;
        OMEGA_FUN(init)(m_parameters.data(), scheme);
    }

    void new_event(const double* data) {OMEGA_FUN(new_event)(data);}
};

struct OmegaWrapperFunctor {
    OmegaWrapper m_omw;

    OmegaWrapperFunctor(OmegaWrapper omw) : m_omw(omw) {}

    double operator()(std::vector<double> lvec, int flv) {
        return m_omw.sqme(lvec.data(), flv);
    } 
};
