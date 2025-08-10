#pragma once
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TH1D.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <Math/MatrixFunctions.h>
#include <ranges>

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


template <unsigned short n_obs>
struct test_struct {
    using Mat_t = ROOT::Math::SMatrix<double, n_obs, n_obs, ROOT::Math::MatRepSym<double, n_obs>>;
    std::vector<Mat_t> m_C;

};

template <unsigned short n_obs>
struct fit_fun {
    // const unsigned short m_n_obs;
    const unsigned short m_n_couplings;
    // pols: LL, LR, RL, RR
    const double m_process_e_pol[4] = {-1., -1., 1., 1.,};
    const double m_process_p_pol[4] = {-1., 1., -1., 1.,};
    // XXX: this could probably be better than just a bunch of vectors...
    // but passing a vector from python is trivial while arrays are annoying
    std::vector<double> m_template_lumi;
    std::vector<std::vector<std::vector<std::vector<double>>>> m_template_param;
    std::vector<std::vector<std::vector<double>>> m_bin_contents_sm;
    // assume all obs having the same binning for now, but can be changed
    // already int in ROOT...
    const int m_n_bins;
    std::vector<double> m_bin_midpoints;
    using Mat_t = ROOT::Math::SMatrix<double, n_obs, n_obs, ROOT::Math::MatRepSym<double, n_obs>>;
    using Vec_t = ROOT::Math::SVector<double, n_obs>;
    using MatAsVec_t = ROOT::Math::SVector<double, n_obs * (n_obs + 1) / 2>;
    std::vector<Mat_t> m_C;

    fit_fun(/*const unsigned short n_obs,*/ const unsigned short n_couplings, const int n_bins, std::vector<double> template_lumi,
        std::vector<std::vector<std::vector<std::vector<double>>>> template_param,
        std::vector<std::vector<std::vector<double>>> bin_contents_sm,
        std::vector<double> bin_midpoints,
        std::vector<MatAsVec_t> C
    ) : m_n_couplings(n_couplings), m_n_bins(n_bins), m_template_lumi(template_lumi),
        m_template_param(template_param),
        m_bin_contents_sm(bin_contents_sm),
        m_bin_midpoints(bin_midpoints)
    {
        std::cout << "start constructor" << std::endl;
        m_C.reserve(4);
        for (const auto& mvec : C) {
            m_C.push_back(Mat_t(mvec, true));
        }
        std::cout << "finished constructor" << std::endl;
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

        // TODO: rewrite to directly use vectors later!
        double diff[n_obs];
        for (unsigned short i = 0; i < n_obs; i++) {
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
                    // XXX: order is unfortunately different here...
                    double bin_content = m_bin_contents_sm[j][i][l] * lumi_scale[j];
                    helicity_contribution += bin_modifier * bin_content * m_bin_midpoints[l];
                }
                temp += helicity_weight[j] * helicity_contribution;
                // TODO add background subtraction (maybe after refactoring all of this into a method?)
            }
            diff[i] = obs[i] - temp;
        }

        // recalculate C_inv from current values of lumi and polarisation
        // might kill the performance in the case that they are fixed...
        Mat_t C_inv;
        for (unsigned short i = 0; i < 4; i++) {
            C_inv += lumi_scale[i] * helicity_weight[i] * m_C[i];
        }
        // TODO: figure out if fast is good enough
        // looks like it but first debug the rest...
        // C_inv.InvertFast();
        C_inv.Invert();

        Vec_t diff_v(diff, n_obs);

        double chi2 = ROOT::Math::Similarity(diff_v, C_inv);
        // std::cout << "chi2: " << chi2 << std::endl;
        // std::cout << "diff_v: " << diff_v << std::endl;
        // std::cout << "C_inv: " << C_inv << std::endl;
        return chi2;
    }
    // takes m_n_obs many observables
    // and 1 (lumi) + 2 (beam pols) + m_n_couplings many parameters
    double operator()(std::vector<double> obs, std::vector<double> pars) {
        return operator()(obs.data(), pars.data());
    }
};

template <unsigned short n_obs, unsigned short n_bins, unsigned short n_couplings>
struct fit_fun2 {
    using Mat_t = ROOT::Math::SMatrix<double, n_obs, n_obs, ROOT::Math::MatRepSym<double, n_obs>>;
    using Vec_t = ROOT::Math::SVector<double, n_obs>;
    using MatAsVec_t = ROOT::Math::SVector<double, n_obs * (n_obs + 1) / 2>;
    // XXX: needs to be scaled to 1 fb^-1
    std::vector<Mat_t> m_C;

    // requires all observables to have the same binning...
    using BinVec_t = ROOT::Math::SVector<double, n_bins>;
    using CouplingVec_t = ROOT::Math::SVector<double, n_couplings>;
    BinVec_t m_binCenters;
    // std::array<std::vector<BinVec_t>, 4> m_templatePar;
    using TemplMat_t = ROOT::Math::SMatrix<double, n_couplings, n_bins>;
    std::array<std::array<TemplMat_t, 4>, n_obs> m_templatePars;

    // XXX: needs to be scaled to 1 fb^-1
    std::array<std::array<BinVec_t, 4>, n_obs> m_signalHists;
    std::vector<std::array<std::array<BinVec_t, 4>, n_obs>> m_backgroundHists;

    // pols: LL, LR, RL, RR
    const double m_process_e_pol[4] = {-1., -1., 1., 1.,};
    const double m_process_p_pol[4] = {-1., 1., -1., 1.,};

    fit_fun2(
        std::vector<std::vector<TH1D*>> signalHists,
        std::vector<double> signalSampleLumi,
        std::vector<MatAsVec_t> C,
        std::vector<std::vector<std::vector<TH1D*>>> templatePars
    ) {
        std::cout << "start constructor" << std::endl;
        // sanity checks of the vectors
        if (signalHists.size() != n_obs) {
            std::cout << "ERROR: provided signalHists.size() == " << signalHists.size() << " but n_obs == " << n_obs << std::endl;
            return;
        }
        if (signalSampleLumi.size() != 4) {
            std::cout << "ERROR: provided signalSampleLumi.size() == " << signalSampleLumi.size() << "but should be 4" << std::endl;
            return;
        }
        if (templatePars.size() != n_obs) {
            std::cout << "ERROR: provided templatePars.size() == " << templatePars.size() << " but n_obs == " << n_obs << std::endl;
            return;
        }
        if (templatePars[0].size() != 4) {
            std::cout << "ERROR: provided templatePars[0].size() == " << templatePars[0].size() << " but should be 4" << std::endl;
            return;
        }
        if (templatePars[0][0].size() != n_couplings) {
            std::cout << "ERROR: provided templatePars[0][0].size() == " << templatePars[0][0].size() << " but n_couplings == " << n_couplings << std::endl;
            return;
        }

        // take first histo to set bin centers
        TH1D* h = signalHists[0][0];
        for (unsigned short i_bin = 0; i_bin < n_bins; i_bin++) {
            m_binCenters[i_bin] = h->GetBinCenter(i_bin+1);
        }

        for (unsigned short i_obs = 0; i_obs < n_obs; i_obs++) {
            for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
                // convert signal histograms
                TH1D* hist = signalHists[i_obs][j_hel];
                // + 1 to skip underflow bin
                BinVec_t signal_vec(hist->GetArray() + 1, n_bins);
                // divide out sample luminosity!
                // FIXME: this should really be done externally!!!
                signal_vec /= signalSampleLumi[j_hel];
                m_signalHists[i_obs][j_hel] = signal_vec;


                // convert templatePars
                TemplMat_t& mat = m_templatePars[i_obs][j_hel];
                // TemplMat_t mat{};
                for (unsigned short k_cpl = 0; k_cpl < n_couplings; k_cpl++) {
                    TH1D* hist = templatePars[i_obs][j_hel][k_cpl];
                    BinVec_t vec(hist->GetArray() + 1, n_bins);
                    mat.Place_in_row(vec, k_cpl, 0);
                    // for (unsigned short l = 0; l < n_bins; l++) {
                    //     mat[k_cpl][l] = vec[l];
                    // }
                }
                // m_templatePars[i_obs][j_hel] = mat;
            }
        }
        m_C.reserve(4);
        for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
            // Mat_t mat(C[j_hel], true);
            Mat_t mat(C[j_hel], false);
            mat /= signalSampleLumi[j_hel];
            m_C.push_back(mat);
        }
        std::cout << "finished constructor" << std::endl;
    }

    fit_fun2(
        std::vector<std::vector<TH1D*>> signalHists,
        std::vector<double> signalSampleLumi,
        std::vector<MatAsVec_t> C,
        std::vector<std::vector<std::vector<TH1D*>>> templatePars,
        std::vector<std::vector<std::vector<TH1D*>>> backgroundHists
    ) : fit_fun2(signalHists, signalSampleLumi, C, templatePars)
    {
        // parse only the background histograms here and do the rest in the original constructor
        for (const auto& bkg : backgroundHists) {
            // check that we have an entry per obs
            if (bkg.size() != n_obs) {
                std::cout << "ERROR: provided bkg.size() == " << bkg.size() << " but n_obs == " << n_obs << std::endl;
                return;
            }
            std::array<std::array<BinVec_t, 4>, n_obs> tmp;
            for (unsigned short i_obs = 0; i_obs < n_obs; i_obs++) {
                // check that we have an entry per helicity
                const auto& obs = bkg[i_obs];
                if (obs.size() != 4) {
                    std::cout << "ERROR: provided obs.size() == " << obs.size() << " but should be 4" << std::endl;
                    return;
                }
                for (unsigned short j_hel; j_hel < 4; j_hel++) {
                    // convert histograms
                    TH1D* hist = bkg[i_obs][j_hel];
                    // + 1 to skip underflow bin
                    BinVec_t bkg_vec(hist->GetArray() + 1, n_bins);
                    tmp[i_obs][j_hel] = bkg_vec;
                }
            }
            m_backgroundHists.push_back(std::move(tmp));
        }
    }


    Vec_t getExpectedSignal(double lumi, std::array<double, 4> helicity_weights,
                            CouplingVec_t couplings, double norm_factor)
        {
            Vec_t res_v{};
            // for each obs
            for (unsigned short i_obs = 0; i_obs < n_obs; i_obs++) {
                BinVec_t expected;
                // for each helicity
                for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
                    // calculate relative changes for this obs and helicity
                    BinVec_t relativeBinChanges = couplings * m_templatePars[i_obs][j_hel];
                    relativeBinChanges += 1.;
                    // element wise multiplication
                    BinVec_t absoluteBinContent = relativeBinChanges * m_signalHists[i_obs][j_hel];
                    // scale to correct lumi and pol
                    double weight = lumi * helicity_weights[j_hel];
                    expected += weight * absoluteBinContent;
                }
                // do the "integral" and apply signal norm factor
                // std::cout << "expected: " << expected << " bin_centers: " << m_binCenters << std::endl;
                double foo = ROOT::Math::Dot(expected, m_binCenters);
                res_v[i_obs] = ROOT::Math::Dot(expected, m_binCenters) * norm_factor;
                // std::cout << "foo: " << foo << " norm_factor: " << norm_factor << std::endl;
            }
            return res_v;
        }


    Vec_t getExpectedBackground(double lumi, std::array<double, 4> helicity_weights,
                                std::vector<double> norm_factors)
        {
            Vec_t res_v{};
            if (norm_factors.size() != m_backgroundHists.size()) {
                std::cout << "ERROR incorrect number of background norm factors given: " << norm_factors.size() << " (needed: " << m_backgroundHists.size() << ")" << std::endl;
                return res_v;
            }
            // for each bkg
            // for (const auto& [bkg, norm_factor] : std::views::zip(m_backgroundHists, norm_factors)) {
            for (unsigned short k_bkg = 0; k_bkg < m_backgroundHists.size(); k_bkg++) {
                const auto& bkg = m_backgroundHists[k_bkg];
                double norm_factor = norm_factors[k_bkg];
                // for each obs
                for (unsigned short i_obs = 0; i_obs < n_obs; i_obs++) {
                    BinVec_t expected;
                    // for each helicity
                    for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
                        // scale to correct lumi and pol
                        double weight = lumi * helicity_weights[j_hel];
                        expected += weight * bkg[i_obs][j_hel];
                    }
                    // do the "integral" and apply norm factor
                    double foo = ROOT::Math::Dot(expected, m_binCenters);
                    res_v[i_obs] += ROOT::Math::Dot(expected, m_binCenters) * norm_factor;
                }
            }
            return res_v;
        }

    Vec_t getExpected(double lumi, std::array<double, 4> helicity_weights,
                      CouplingVec_t couplings, std::vector<double>norm_factors)
        {
            Vec_t signal = this->getExpectedSignal(lumi, helicity_weights, couplings, norm_factors[0]);

            if (m_backgroundHists.size() > 0) {
                Vec_t background = this->getExpectedBackground(lumi, helicity_weights, norm_factors);
                signal += background;
            }

            return signal;
        }

    double operator()(double* obs, double* pars) {
        // TODO: all of this should probably be const...
        double lumi = pars[0];
        double beam_e_pol = pars[1];
        double beam_p_pol = pars[2];
        // implicitly ties n_obs to n_couplings, which makes sense for OO
        // but I could also fit more couplings than obs if I want...
        CouplingVec_t couplings(pars + 3, n_obs);
        // TODO: only one so far, add more for backgrounds
        std::vector<double> norm_factors(pars + 3 + n_obs, pars + 3 + n_obs + 1);

        // TODO: store last beam pols and cache weights and matrix if fit is too slow
        std::array<double, 4> helicity_weight;
        for (unsigned short i = 0; i < 4; i++) {
            helicity_weight[i] = 0.25 * (1.0 + beam_e_pol * m_process_e_pol[i]) * (1.0 + beam_p_pol * m_process_p_pol[i]);
        }
        // recalculate C_inv from current values of lumi and polarisation
        // might kill the performance in the case that they are fixed...
        Mat_t C_inv;
        for (unsigned short i = 0; i < 4; i++) {
            C_inv += helicity_weight[i] * m_C[i];
        }
        C_inv *= lumi;
        // std::cout << "C: " << C_inv << std::endl;
        // TODO: figure out if fast is good enough
        // looks like it but first debug the rest...
        C_inv.InvertFast();
        // C_inv.Invert();

        Vec_t obs_v(obs, n_obs);
        Vec_t exp_v = this->getExpected(lumi, helicity_weight, couplings, norm_factors);

        Vec_t diff_v = obs_v - exp_v;

        double chi2 = ROOT::Math::Similarity(diff_v, C_inv);
        // std::cout << "chi2: " << chi2 << std::endl;
        // std::cout << "obs_v: " << obs_v << std::endl;
        // std::cout << "exp_v: " << exp_v << std::endl;
        // std::cout << "diff_v: " << diff_v << std::endl;
        // std::cout << "C_inv: " << C_inv << std::endl;
        return chi2;
    }

    // takes m_n_obs many observables
    // and 1 (lumi) + 2 (beam pols) + m_n_couplings many parameters
    // needed for ROOT::Math::IBaseFunctionMultiDim
    double operator()(double* all) {
        double* obs = all;
        double* pars = all + n_obs;
        return operator()(obs, pars);
    }

    // takes m_n_obs many observables
    // and 1 (lumi) + 2 (beam pols) + m_n_couplings many parameters
    // easier to call with cppyy for debugging
    double operator()(std::vector<double> obs, std::vector<double> pars) {
        return operator()(obs.data(), pars.data());
    }
};

template <unsigned short n_bins, unsigned short n_couplings>
struct expect_fun {
    // requires all observables to have the same binning...
    using BinVec_t = ROOT::Math::SVector<double, n_bins>;
    using CouplingVec_t = ROOT::Math::SVector<double, n_couplings>;
    BinVec_t m_binCenters;
    // std::array<std::vector<BinVec_t>, 4> m_templatePar;
    using TemplMat_t = ROOT::Math::SMatrix<double, n_couplings, n_bins>;
    std::array<TemplMat_t, 4> m_templatePars;

    // XXX: needs to be scaled to 1 fb^-1
    std::array<BinVec_t, 4> m_signalHists;
    std::vector<std::array<BinVec_t, 4>> m_backgroundHists;

    // pols: LL, LR, RL, RR
    const double m_process_e_pol[4] = {-1., -1., 1., 1.,};
    const double m_process_p_pol[4] = {-1., 1., -1., 1.,};

    // const double m_lumi;
    // const double m_e_pol;
    // const double m_p_pol;

    expect_fun(
        std::vector<TH1D*> signalHists,
        std::vector<std::vector<TH1D*>> templatePars
        /*double lumi, double e_pol, double p_pol*/
    ) /*: m_lumi(lumi), m_e_pol(e_pol), m_p_pol(p_pol)*/
    {
        std::cout << "start constructor" << std::endl;
        // sanity checks of the vectors
        if (templatePars.size() != 4) {
            std::cout << "ERROR: provided templatePars.size() == " << templatePars.size() << " but should be 4" << std::endl;
            return;
        }
        if (templatePars[0].size() != n_couplings) {
            std::cout << "ERROR: provided templatePars[0].size() == " << templatePars[0].size() << " but n_couplings == " << n_couplings << std::endl;
            return;
        }

        // take first histo to set bin centers
        TH1D* h = signalHists[0];
        for (unsigned short i_bin = 0; i_bin < n_bins; i_bin++) {
            m_binCenters[i_bin] = h->GetBinCenter(i_bin+1);
        }

        for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
            // convert signal histograms
            TH1D* hist = signalHists[j_hel];
            // + 1 to skip underflow bin
            BinVec_t signal_vec(hist->GetArray() + 1, n_bins);
            m_signalHists[j_hel] = signal_vec;

            // convert templatePars
            TemplMat_t& mat = m_templatePars[j_hel];
            for (unsigned short k_cpl = 0; k_cpl < n_couplings; k_cpl++) {
                TH1D* hist = templatePars[j_hel][k_cpl];
                BinVec_t vec(hist->GetArray() + 1, n_bins);
                mat.Place_in_row(vec, k_cpl, 0);
            }
        }
        std::cout << "finished constructor" << std::endl;
    }

    expect_fun(
        std::vector<TH1D*> signalHists,
        std::vector<std::vector<TH1D*>> templatePars,
        /*double lumi, double e_pol, double p_pol,*/
        std::vector<std::vector<TH1D*>> backgroundHists
    ) : expect_fun(signalHists, templatePars/*, lumi, e_pol, p_pol*/)
    {
        // parse only the background histograms here and do the rest in the original constructor
        for (const auto& bkg : backgroundHists) {
            std::array<BinVec_t, 4> tmp;
            // check that we have an entry per helicity
            if (bkg.size() != 4) {
                std::cout << "ERROR: provided bkg.size() == " << bkg.size() << " but should be 4" << std::endl;
                return;
            }
            for (unsigned short j_hel; j_hel < 4; j_hel++) {
                // convert histograms
                TH1D* hist = bkg[j_hel];
                // + 1 to skip underflow bin
                BinVec_t bkg_vec(hist->GetArray() + 1, n_bins);
                tmp[j_hel] = bkg_vec;
            }
        m_backgroundHists.push_back(std::move(tmp));
        }
    }


    double getExpectedSignal(double lumi, std::array<double, 4> helicity_weights,
                            CouplingVec_t couplings, double norm_factor)
        {
            BinVec_t expected;
            // for each helicity
            for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
                // calculate relative changes for this obs and helicity
                BinVec_t relativeBinChanges = couplings * m_templatePars[j_hel];
                relativeBinChanges += 1.;
                // element wise multiplication
                BinVec_t absoluteBinContent = relativeBinChanges * m_signalHists[j_hel];
                // scale to correct lumi and pol
                double weight = lumi * helicity_weights[j_hel];
                expected += weight * absoluteBinContent;
            }
            double foo = ROOT::Math::Dot(expected, m_binCenters);
            return ROOT::Math::Dot(expected, m_binCenters) * norm_factor;
        }


    double getExpectedBackground(double lumi, std::array<double, 4> helicity_weights,
                                std::vector<double> norm_factors)
        {
            double res = 0.;
            if (norm_factors.size() != m_backgroundHists.size()) {
                std::cout << "ERROR incorrect number of background norm factors given: " << norm_factors.size() << " (needed: " << m_backgroundHists.size() << ")" << std::endl;
                return res;
            }
            // for each bkg
            // for (const auto& [bkg, norm_factor] : std::views::zip(m_backgroundHists, norm_factors)) {
            for (unsigned short k_bkg = 0; k_bkg < m_backgroundHists.size(); k_bkg++) {
                const auto& bkg = m_backgroundHists[k_bkg];
                double norm_factor = norm_factors[k_bkg];
                BinVec_t expected;
                // for each helicity
                for (unsigned short j_hel = 0; j_hel < 4; j_hel++) {
                    // scale to correct lumi and pol
                    double weight = lumi * helicity_weights[j_hel];
                    expected += weight * bkg[j_hel];
                }
                // do the "integral" and apply norm factor
                double foo = ROOT::Math::Dot(expected, m_binCenters);
                res += ROOT::Math::Dot(expected, m_binCenters) * norm_factor;
            }
            return res;
        }

    double getExpected(double lumi, std::array<double, 4> helicity_weights,
                      CouplingVec_t couplings, std::vector<double>norm_factors)
        {
            double signal = this->getExpectedSignal(lumi, helicity_weights, couplings, norm_factors[0]);

            if (m_backgroundHists.size() > 0) {
                double background = this->getExpectedBackground(lumi, helicity_weights, norm_factors);
                signal += background;
            }

            return signal;
        }

    double operator()(double* pars) {
        // TODO: all of this should probably be const...
        double lumi = pars[0];
        double beam_e_pol = pars[1];
        double beam_p_pol = pars[2];
        CouplingVec_t couplings(pars + 3, n_couplings);
        // TODO: only one so far, add more for backgrounds
        std::vector<double> norm_factors(pars + 3 + n_couplings, pars + 3 + n_couplings + 1);

        std::array<double, 4> helicity_weight;
        for (unsigned short i = 0; i < 4; i++) {
            helicity_weight[i] = 0.25 * (1.0 + beam_e_pol * m_process_e_pol[i]) * (1.0 + beam_p_pol * m_process_p_pol[i]);
        }

        return this->getExpected(lumi, helicity_weight, couplings, norm_factors);
    }

    // easier to call with cppyy for debugging
    double operator()(std::vector<double> pars) {
        return operator()(pars.data());
    }
};
