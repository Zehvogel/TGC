#pragma once
#include <root/TMatrixDSym.h>
#include <root/TVectorD.h>
#include <root/TH1D.h>
#include <root/Math/SMatrix.h>
#include <root/Math/SVector.h>
#include <root/Math/MatrixFunctions.h>

// struct fit_fun {
//     unsigned short m_n_obs;
//     TMatrixDSym m_C_inv_t;

//     fit_fun(unsigned short n_obs, TMatrixDSym& C) : m_n_obs(n_obs), m_C_inv_t(TMatrixDSym::EMatrixCreatorsOp1::kTransposed, TMatrixDSym(TMatrixDSym::EMatrixCreatorsOp1::kInverted, C)) {}

//     double operator()(double* obs, double* pars) {
//         double diff[m_n_obs];
//         for (unsigned short i = 0; i < m_n_obs; i++) {
//             diff[i] = obs[i] - template_for_pars(pars);
//         }
//         TVectorD v(m_n_obs, diff);
//         TVectorD C_inv_t_v = m_C_inv_t * v;
//         double res = 0;
//         for (unsigned short i = 0; i < m_n_obs; i++) {
//             res += C_inv_t_v[i] * v[i];
//         }
//         return res;
//     }
// };

// struct hist_from_template {
//     // FIXME: maybe just get rid of histograms and use vectors instead!
//     // I need to add up all the stuff in the end anyway so I could also just do it directly in another loop!
//     std::vector<TH1D&> m_hists;
//     // store SM hist in first entry
//     hist_from_template(std::vector<TH1D&> hists) : m_hists(hists) {}

//     TH1D operator()(std::vector<double> parameters) {
//         TH1D* res = dynamic_cast<TH1D*>(m_hists[0].Clone());
//         for (size_t i = 0; i < parameters.size(); i++) {
//             *res += parameters[i] * m_hists[i+1];
//         }
//         return *res;
//     }
// };

struct fit_fun {
    const unsigned short m_n_obs;
    const unsigned short m_n_couplings;
    // already int in ROOT...
    const int m_n_bins;
    // pols: LL, LR, RL, RR
    const double m_process_e_pol[4] = {-1., -1., 1., 1.,};
    const double m_process_p_pol[4] = {-1., 1., -1., 1.,};
    // XXX: this could probably be better than just a bunch of vectors...
    std::vector<const double> m_template_lumi;
    std::vector<std::vector<std::vector<std::vector<const double>>>> m_template_param;
    std::vector<std::vector<std::vector<const double>>> m_bin_contents_sm;
    std::vector<const double> m_bin_midpoints;

    fit_fun(const unsigned short n_obs, const unsigned short n_couplings, const int n_bins, std::vector<const double> template_lumi,
        std::vector<std::vector<std::vector<std::vector<const double>>>> template_param,
        std::vector<std::vector<std::vector<const double>>> bin_contents_sm,
        std::vector<const double> bin_midpoints
    ) : m_n_obs(n_obs), m_n_couplings(n_couplings), m_n_bins(n_bins), m_template_lumi(template_lumi)
    {
        // empty?
    }

    // takes m_n_obs many observables
    // and 1 (lumi) + 2 (beam pols) + m_n_couplings many parameters
    double operator()(double* obs, double* pars) {
        // TODO: all of this should probably be const...
        double lumi = pars[0];
        double beam_e_pol = pars[1];
        double beam_p_pol = pars[2];
        double* couplings = pars + 3;

        double lumi_scale[4] = {0.};
        double helicity_weight[4] = {0.};
        for (unsigned short i = 0; i < 4; i++) {
            lumi_scale[i] = lumi / m_template_lumi[i];
            helicity_weight[i] = 0.25 * (1.0 + beam_e_pol * m_process_e_pol[i]) * (1.0 + beam_p_pol * m_process_p_pol[i]);
        }

        double diff[m_n_obs];
        for (unsigned short i = 0; i < m_n_obs; i++) {
            double temp = 0;
            // loop over helicity combinations
            for (int j = 0; j < 4; j++) {
                // loop over bins
                double helicity_contribution = 0;
                for (int l = 0; l < m_n_bins; l++) {
                    double bin_modifier = 1.;
                    // loop over couplings
                    for (unsigned short k = 0; k < m_n_couplings; k++) {
                        // TODO add mixed and quadratic terms
                        bin_modifier += couplings[k] * m_template_param[i][j][k][l];
                    }
                    double bin_content = m_bin_contents_sm[i][j][l] * lumi_scale[j];
                    helicity_contribution += bin_modifier * bin_content * m_bin_midpoints[l];
                }
                temp += helicity_weight[j] * helicity_contribution;
                // TODO add background subtraction (maybe after refactoring all of this into a method?)
            }
            diff[i] = obs[i] - temp;
        }

        // recalculate C_inv from current values of lumi and polarisation
        // might kill the performance a bit in the case that they are fixed...
        unsigned int n = m_n_obs * (m_n_obs + 1) / 2;
        double C_inv[n];
        for (unsigned int i = 0; i < n; i++) {
            double c = 0;
            for (unsigned short j = 0; j < 4; j++) {
                // inverse lumi scale as we want C~_inv with C~ = lumi*C
                // FIXME: is inversion of a matrix a linear operation?? I don't think so!
                // So I would need to assemble a matrix C from the helicities and then invert it in every call???????
                // that will certainly slow things down...
                // TODO: change it all to use SMatrix and SVector, I will not implement cramers rule for n_obs = 6 myself!
                c = (1 / lumi_scale[j]) * helicity_weight[j] * m_C_inv[j][i];
            }
            C_inv[i] = c;
        }

        double chi2 = 0.;
        unsigned short k = 0;
        for (unsigned short i = 0; i < m_n_obs; i++) {
            // C_inv is symmetric so we can save a bit in the inner loop
            for (unsigned short j = 0; j <= i; j++, k++) {
                double tmp = diff[i] * C_inv[k] * diff[j];
                if (i != j) {
                    // we need to add the contribution from the off-diagonals twice.
                    tmp *= 2.;
                }
                chi2 += tmp;
            }
        }
    }
};
