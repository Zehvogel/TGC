from analysis_framework.Analysis import Analysis
import ROOT
import numpy as np
from OO.whizard.model_parser import ModelParser
import subprocess


def make_lvec_E(column: str, idx: str):
    return f"""
    ROOT::Math::PxPyPzEVector(
        {column}.momentum.x[{idx}],
        {column}.momentum.y[{idx}],
        {column}.momentum.z[{idx}],
        {column}.energy[{idx}]
    )
    """


def make_lvec_M(column: str, idx: str):
    return f"""
    ROOT::Math::PxPyPzMVector(
        {column}.momentum.x[{idx}],
        {column}.momentum.y[{idx}],
        {column}.momentum.z[{idx}],
        {column}.mass[{idx}]
    )
    """


class WWAnalysis(Analysis):

    truth_defined: bool
    truth_categories: list[str]
    _omega_wrappers = {}
    _mc_indices = {}
    _signal_categories: list[str]

    def __init__(self, dataset):
        self.truth_defined = False
        self.truth_categories = []
        super().__init__(dataset)


    def define_reco_objects(self, x_angle: float):
        # select isolated lepton and the two jets and build Ws and neutrino
        self.Define("iso_lep_idx", "IsolatedElectrons_objIdx.index[0]")
        self.Define("iso_lep_charge", "PandoraPFOs.charge[iso_lep_idx]")
        # TODO: figure out: are isolated electrons also always Pandora electrons? :( (taking a look into the IsolatedLeptonTagging shows that this is not required!)
        self.Define("iso_lep_lvec", "ROOT::Math::PxPyPzEVector(PandoraPFOs.momentum.x[iso_lep_idx], PandoraPFOs.momentum.y[iso_lep_idx], PandoraPFOs.momentum.z[iso_lep_idx], PandoraPFOs.energy[iso_lep_idx])")
        self.Define("jet1_lvec", "R2jet_lvecs[0]")
        self.Define("jet2_lvec", "R2jet_lvecs[1]")
        self.Define("sqrt_s_E", "params_Energy")
        # XXX: ignoring the electron mass here
        self.Define("sqrt_s_lvec", f"ROOT::Math::PxPyPzEVector(sqrt_s_E * sin({x_angle}/2.), 0, 0, sqrt_s_E)")
        self.Define("hadronic_W_lvec", "jet1_lvec + jet2_lvec")
        # TODO: overlay contamination etc are still in this...
        self.Define("leptonic_W_lvec", "sqrt_s_lvec - hadronic_W_lvec")
        self.Define("nu_lvec", "leptonic_W_lvec - iso_lep_lvec")


    def remove_x_angle(self, x_angle: float):
        # now we would want to remove the crossing angle
        #small experiment using an actual dd4hep crossing angle boosted e+e- pair:
        # e_lvec = ROOT.Math.PxPyPzMVector(+8.750143e-01, 0., +1.250000e+02, +5.109968e-04)
        # p_lvec = ROOT.Math.PxPyPzMVector(+8.750143e-01, 0., -1.250000e+02, +5.109968e-04)
        # s_lvec = e_lvec + p_lvec
        # print(s_lvec)
        # beta = s_lvec.BoostToCM()
        # print(beta.x(), beta.y(), beta.z())
        # boost = ROOT.Math.BoostX(beta.x())
        # print(boost(s_lvec))
        # print(boost(e_lvec))
        # print(ROOT.Math.sin(0.007))
        # boost2 = ROOT.Math.BoostX(-ROOT.Math.sin(0.007))
        # print(boost2(s_lvec))
        # print(boost2(e_lvec))
        ROOT.gInterpreter.Declare(f"ROOT::Math::BoostX unboost_xangle(-std::sin({x_angle/2}));")
        lvec_list = [
            "iso_lep_lvec", "jet1_lvec", "jet2_lvec",
            "hadronic_W_lvec", "leptonic_W_lvec", "nu_lvec"
            ]
        for lvec in lvec_list:
            self.Define(f"ub_{lvec}", f"unboost_xangle({lvec})")
        if self.truth_defined:
            truth_lvec_list = [
                "true_lep_lvec", "true_nu_lvec", "true_quark1_lvec", "true_quark2_lvec",
                "true_leptonic_W_lvec", "true_hadronic_W_lvec"
                ]
            for lvec in truth_lvec_list:
                self.define_only_on(self.truth_categories, f"ub_{lvec}", f"unboost_xangle({lvec})")


    def define_hagiwara_angles(self):
        # now we want to calculate the angles as in the Hagiwara parametrisation
        # XXX: there is this one decision to make, do we use both measured daughter values or just one and flip it? :(
        # TODO: is this well defined for the e+nuqq case?
        # I need to be a bit careful, there can be some replacements that are only valid in zero-width CC03
        # It all depends on the systems of the measurement actually being the rest-frame of the parent i.e. that the decay are back-to-back
        # Exactly in this condition ang(W-, e-) == ang(W+, e+) == Theta_W
        # In the W- rest-frame f has theta and phi and fbar has -theta and phi, what are these angles actually
        # theta ang(f*, z') with f* in a system rotated such that z' direction of W- in the CMS and y' = P_e- x P_W-
        # boosted into the rest-frame of the parent of f
        # for now we just act as if we are in the double pole case as in Hagiwara and proceed.
        # TODO: guaranteed to not get 0 charge?
        self.Define("Wm_lvec", "iso_lep_charge < 0. ? ub_leptonic_W_lvec : ub_hadronic_W_lvec")
        self.Define("Wp_lvec", "iso_lep_charge > 0. ? ub_leptonic_W_lvec : ub_hadronic_W_lvec")
        # they literally did not put a cheaper cosTheta accessor into Genvectors...
        self.Define("Wm_cosTheta", "cos(Wm_lvec.Theta())")

        ROOT.gInterpreter.Declare("#include <WWTools.h>")
        self.Define("iso_lep_star_lvec", "WWTools::starVectorHagiwara(ub_leptonic_W_lvec, ub_iso_lep_lvec, {0, 0, 1}, Wm_lvec)")
        self.Define("nu_star_lvec", "WWTools::starVectorHagiwara(ub_leptonic_W_lvec, ub_nu_lvec, {0, 0, 1}, Wm_lvec)")

        self.Define("jet1_star_lvec", "WWTools::starVectorHagiwara(ub_hadronic_W_lvec, ub_jet1_lvec, {0, 0, 1}, Wm_lvec)")
        self.Define("jet2_star_lvec", "WWTools::starVectorHagiwara(ub_hadronic_W_lvec, ub_jet2_lvec, {0, 0, 1}, Wm_lvec)")

        lvec_list = ["iso_lep", "nu", "jet1", "jet2"]
        for lvec in lvec_list:
            self.Define(f"{lvec}_co", f"cos({lvec}_star_lvec.Theta())")
            self.Define(f"{lvec}_ph", f"{lvec}_star_lvec.Phi()")

        if self.truth_defined:
            self.define_only_on(self.truth_categories, "true_Wm_lvec", "true_iso_lep_charge < 0. ? ub_true_leptonic_W_lvec : ub_true_hadronic_W_lvec")
            self.define_only_on(self.truth_categories, "true_Wp_lvec", "true_iso_lep_charge > 0. ? ub_true_leptonic_W_lvec : ub_true_hadronic_W_lvec")
            self.define_only_on(self.truth_categories, "true_Wm_cosTheta", "cos(true_Wm_lvec.Theta())")

            self.define_only_on(self.truth_categories, "true_iso_lep_star_lvec", "WWTools::starVectorHagiwara(ub_true_leptonic_W_lvec, ub_true_lep_lvec, {0, 0, 1}, true_Wm_lvec)")
            self.define_only_on(self.truth_categories, "true_nu_star_lvec", "WWTools::starVectorHagiwara(ub_true_leptonic_W_lvec, ub_true_nu_lvec, {0, 0, 1}, true_Wm_lvec)")
            self.define_only_on(self.truth_categories, "true_jet1_star_lvec", "WWTools::starVectorHagiwara(ub_true_hadronic_W_lvec, ub_true_quark1_lvec, {0, 0, 1}, true_Wm_lvec)")
            self.define_only_on(self.truth_categories, "true_jet2_star_lvec", "WWTools::starVectorHagiwara(ub_true_hadronic_W_lvec, ub_true_quark2_lvec, {0, 0, 1}, true_Wm_lvec)")
            truth_lvec_list = ["true_iso_lep", "true_nu", "true_jet1", "true_jet2"]
            for lvec in truth_lvec_list:
                self.define_only_on(self.truth_categories, f"{lvec}_co", f"cos({lvec}_star_lvec.Theta())")
                self.define_only_on(self.truth_categories, f"{lvec}_ph", f"{lvec}_star_lvec.Phi()")



    def define_OO(self):
        self.Define("co", "Wm_cosTheta")
        ROOT.gInterpreter.Declare("""
        double proba_fold_Wm(double co, double co11, double co12, double co2, double ph11, double ph12, double ph2, int cpl, int order)
        {
                                  double res = 0.;
                                  res += OOTools::proba(co, co11, co2, ph11, ph2, cpl, false, order);
                                  res += OOTools::proba(co, co12, co2, ph12, ph2, cpl, false, order);
                                  return res;
        }
        double proba_fold_Wp(double co, double co1, double co21, double co22, double ph1, double ph21, double ph22, int cpl, int order)
        {
                                  double res = 0.;
                                  res += OOTools::proba(co, co1, co21, ph1, ph21, cpl, false, order);
                                  res += OOTools::proba(co, co1, co22, ph1, ph22, cpl, false, order);
                                  return res;
        }
        """)

        self.Define("S_0",   "iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 0, 0) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 0, 0)")
        self.Define("S_1_1", "iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 1, 1) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 1, 1)")
        self.Define("S_1_2", "iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 2, 1) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 2, 1)")
        self.Define("S_1_3", "iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 3, 1) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 3, 1)")

        self.Define("O_1", "S_1_1 / S_0")
        self.Define("O_2", "S_1_2 / S_0")
        self.Define("O_3", "S_1_3 / S_0")

        # for debug
        self.Define("co1", "iso_lep_charge < 0. ? iso_lep_co : -iso_lep_co")
        self.Define("ph1", "iso_lep_charge < 0. ? iso_lep_ph : iso_lep_ph + ROOT::Math::Pi() <= ROOT::Math::Pi() ? iso_lep_ph + ROOT::Math::Pi() : iso_lep_ph - ROOT::Math::Pi()")

        if self.truth_defined:
            self.define_only_on(self.truth_categories, "true_co", "true_Wm_cosTheta")
            self.define_only_on(self.truth_categories, "true_S_0",   "true_iso_lep_charge < 0. ? proba_fold_Wp(true_co, true_iso_lep_co, true_jet1_co, true_jet2_co, true_iso_lep_ph, true_jet1_ph, true_jet2_ph, 0, 0) : proba_fold_Wm(true_co, true_jet1_co, true_jet2_co, true_nu_co, true_jet1_ph, true_jet2_ph, true_nu_ph, 0, 0)")
            for i in range(1, 4):
                self.define_only_on(self.truth_categories, f"true_S_1_{i}", f"true_iso_lep_charge < 0. ? proba_fold_Wp(true_co, true_iso_lep_co, true_jet1_co, true_jet2_co, true_iso_lep_ph, true_jet1_ph, true_jet2_ph, {i}, 1) : proba_fold_Wm(true_co, true_jet1_co, true_jet2_co, true_nu_co, true_jet1_ph, true_jet2_ph, true_nu_ph, {i}, 1)")
                self.define_only_on(self.truth_categories, f"true_O_{i}", f"true_S_1_{i} / true_S_0")
            # for debug
            self.define_only_on(self.truth_categories, "true_co1", "true_iso_lep_charge < 0. ? true_iso_lep_co : -true_iso_lep_co")
            self.define_only_on(self.truth_categories, "true_ph1", "true_iso_lep_charge < 0. ? true_iso_lep_ph : true_iso_lep_ph + ROOT::Math::Pi() <= ROOT::Math::Pi() ? true_iso_lep_ph + ROOT::Math::Pi() : true_iso_lep_ph - ROOT::Math::Pi()")


    def define_truth_objects(self, categories: list[str]):
        # take first genstat 1 e and nu and first two gen stat 2 pdg below 6
        def first_stable(abs_pdg):
            return f"""
            auto& genStat = MCParticlesSkimmed.generatorStatus;
            auto& pdg = MCParticlesSkimmed.PDG;
            auto mask = genStat == 1 && abs(pdg) == {abs_pdg};
            // abuse ArgMax to get the first set position
            auto idx = ArgMax(mask);
            return idx;
            """
        def first_two_unstable_below(abs_pdg):
            return f"""
            auto& genStat = MCParticlesSkimmed.generatorStatus;
            auto& pdg = MCParticlesSkimmed.PDG;
            auto mask = genStat == 2 && abs(pdg) < {abs_pdg};
            // abuse ArgMax to get the first set position
            auto idx = ArgMax(mask);
            // Drop the first index and use ArgMax again
            auto idx2 = ArgMax(Drop(mask, {{idx}}));
            // increment by one to compensate for removal of the first
            idx2++;
            return ROOT::VecOps::RVec({{idx, idx2}});
            """
        # FIXME: this selection will return over optimistic results as it will compare the MC electron after FSR with the reconstructed one
        # making the reconstructed one appear less wrong than it is, more correct for the purpose of OO would be to take it directly after the ME calc
        if not categories:
            self.Define("true_lep_idx", first_stable(11))
            self.Define("true_nu_idx", first_stable(12))
            self.Define("true_quarks_idcs", first_two_unstable_below(6))
            self.Define("true_quark1_idx", "true_quarks_idcs[0]")
            self.Define("true_quark2_idx", "true_quarks_idcs[1]")
            self.Define("true_lep_lvec", make_lvec_M("MCParticlesSkimmed", "true_lep_idx"))
            self.Define("true_nu_lvec", make_lvec_M("MCParticlesSkimmed", "true_nu_idx"))
            self.Define("true_quark1_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark1_idx"))
            self.Define("true_quark2_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark2_idx"))
            self.Define("true_leptonic_W_lvec", "true_lep_lvec + true_nu_lvec")
            self.Define("true_hadronic_W_lvec", "true_quark1_lvec + true_quark2_lvec")
            self.Define("true_iso_lep_charge", "MCParticlesSkimmed.PDG[true_lep_idx] > 0. ? -1. : 1.")
        else:
            self.define_only_on(categories, "true_lep_idx", first_stable(11))
            self.define_only_on(categories, "true_nu_idx", first_stable(12))
            self.define_only_on(categories, "true_quarks_idcs", first_two_unstable_below(6))
            self.define_only_on(categories, "true_quark1_idx", "true_quarks_idcs[0]")
            self.define_only_on(categories, "true_quark2_idx", "true_quarks_idcs[1]")
            self.define_only_on(categories, "true_lep_lvec", make_lvec_M("MCParticlesSkimmed", "true_lep_idx"))
            self.define_only_on(categories, "true_nu_lvec", make_lvec_M("MCParticlesSkimmed", "true_nu_idx"))
            self.define_only_on(categories, "true_quark1_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark1_idx"))
            self.define_only_on(categories, "true_quark2_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark2_idx"))
            self.define_only_on(categories, "true_leptonic_W_lvec", "true_lep_lvec + true_nu_lvec")
            self.define_only_on(categories, "true_hadronic_W_lvec", "true_quark1_lvec + true_quark2_lvec")
            self.define_only_on(categories, "true_iso_lep_charge", "MCParticlesSkimmed.PDG[true_lep_idx] > 0. ? -1. : 1.")

        self.truth_defined = True
        self.truth_categories = categories


    def book_OO_matrix(self):
        for i in range(1, 4):
            for j in range(1, 4):
                self.Define(f"c_{i}{j}", f"O_{i} * O_{j}")
                self.book_sum(f"c_{i}{j}", f"c_{i}{j}")
                if self.truth_defined:
                    self.define_only_on(self.truth_categories, f"true_c_{i}{j}", f"true_O_{i} * true_O_{j}")
                    self.book_sum(f"true_c_{i}{j}", f"true_c_{i}{j}", categories=self.truth_categories)



    def get_OO_matrix_normalized(self, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0):
        values = [self.get_mean(f"c_{i}{j}", int_lumi, e_pol, p_pol) for i in range(1,4) for j in range(1,4)]
        c = np.asarray(values).reshape((3, 3))
        return c


    def get_true_OO_matrix_normalized(self, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0):
        values = [self.get_mean(f"true_c_{i}{j}", int_lumi, e_pol, p_pol) for i in range(1,4) for j in range(1,4)]
        c = np.asarray(values).reshape((3, 3))
        return c


    def initialise_omega_wrappers(self, configurations: dict[str,dict[str, float]]):

        whizard_prefix = subprocess.run(['whizard-config', '--prefix'], capture_output=True, encoding='ascii').stdout.strip()
        whizard_libs = f"{whizard_prefix}/lib/"
        # print(whizard_libs)
        ROOT.gSystem.AddDynamicPath(whizard_libs)
        ROOT.gSystem.Load("libwhizard.so")
        ROOT.gSystem.Load("libwhizard_main.so")
        ROOT.gSystem.Load("libomega.so")
        ROOT.gSystem.Load("OO/whizard/cc20_ac_inclusive/.libs/default_lib.so")
        ROOT.gInterpreter.Declare("#include \"OO/whizard/OmegaWrapper.h\"")

        model_parser = ModelParser("OO/whizard/SM_ac.mdl")
        # add derivation of lz and kz according to lep parametrisation
        model_parser.add_derived_parameter("lz", "la")
        model_parser.add_derived_parameter("kz", "1.0 - (ka - 1.0) * sw**2/cw**2 + (g1z - 1.0)")
        self._omega_wrappers["nominal"] = ROOT.OmegaWrapper(model_parser.get_parameters_list())

        for name, pars in configurations.items():
            model_parser.set_parameters(pars)
            self._omega_wrappers[name] = ROOT.OmegaWrapper(model_parser.get_parameters_list())


    def set_mc_indices(self, indices: dict[str, int]):
        self._mc_indices = indices


    def set_signal_categories(self, categories: list[str]):
        self._signal_categories = categories


    def calc_reco_sqme(self):
        self.Define("reco_ME_flv", "iso_lep_charge > 0 ? 1 : 2")
        self.Define("reco_ME_momenta_12", """
            std::vector<double>({
                    125., 0., 0., 125.,
                    125., 0., 0., -125.,
                    ub_iso_lep_lvec.E(), ub_iso_lep_lvec.Px(), ub_iso_lep_lvec.Py(), ub_iso_lep_lvec.Pz(),
                    ub_nu_lvec.E(), ub_nu_lvec.Px(), ub_nu_lvec.Py(), ub_nu_lvec.Pz(),
                    ub_jet1_lvec.E(), ub_jet1_lvec.Px(), ub_jet1_lvec.Py(), ub_jet1_lvec.Pz(),
                    ub_jet2_lvec.E(), ub_jet2_lvec.Px(), ub_jet2_lvec.Py(), ub_jet2_lvec.Pz(),
            })
        """)
        self.Define("reco_ME_momenta_21", """
            std::vector<double>({
                    125., 0., 0., 125.,
                    125., 0., 0., -125.,
                    ub_iso_lep_lvec.E(), ub_iso_lep_lvec.Px(), ub_iso_lep_lvec.Py(), ub_iso_lep_lvec.Pz(),
                    ub_nu_lvec.E(), ub_nu_lvec.Px(), ub_nu_lvec.Py(), ub_nu_lvec.Pz(),
                    ub_jet2_lvec.E(), ub_jet2_lvec.Px(), ub_jet2_lvec.Py(), ub_jet2_lvec.Pz(),
                    ub_jet1_lvec.E(), ub_jet1_lvec.Px(), ub_jet1_lvec.Py(), ub_jet1_lvec.Pz(),
            })
        """)
        for name, omw in self._omega_wrappers.items():
            self.Define(f"reco_sqme_12_{name}", omw, ["reco_ME_momenta_12", "reco_ME_flv"])
            self.Define(f"reco_sqme_21_{name}", omw, ["reco_ME_momenta_21", "reco_ME_flv"])


    # FIXME: urgh only do all this for the signal categories
    def book_weights(self):
        self.define_only_on(self._signal_categories, "lep_charge", "MCParticlesSkimmed.charge[true_lep_idx]")
        self.define_only_on(self._signal_categories, "mc_ME_flv", "lep_charge > 0 ? 1 : 2")
        # TODO: optimise
        self.define_only_on(self._signal_categories, "mc_lvec", "Construct<ROOT::Math::PxPyPzMVector>(MCParticlesSkimmed.momentum.x, MCParticlesSkimmed.momentum.y, MCParticlesSkimmed.momentum.z, MCParticlesSkimmed.mass)")
        self.define_only_on(self._signal_categories, "mc_E", "return Map(mc_lvec, [] (const auto& el) {return el.energy();})")
        self.define_only_on(self._signal_categories, "mc_PX", "MCParticlesSkimmed.momentum.x")
        self.define_only_on(self._signal_categories, "mc_PY", "MCParticlesSkimmed.momentum.y")
        self.define_only_on(self._signal_categories, "mc_PZ", "MCParticlesSkimmed.momentum.z")
        beam_e_idx = self._mc_indices["beam_e_ISR"]
        beam_p_idx = self._mc_indices["beam_p_ISR"]
        self.define_only_on(self._signal_categories, "mc_ME_momenta", f"""
                    std::vector<double>({{
                    mc_E[{beam_e_idx}],    mc_PX[{beam_e_idx}],    mc_PY[{beam_e_idx}],    mc_PZ[{beam_e_idx}],
                    mc_E[{beam_p_idx}],    mc_PX[{beam_p_idx}],    mc_PY[{beam_p_idx}],    mc_PZ[{beam_p_idx}],
                    mc_E[true_lep_idx],    mc_PX[true_lep_idx],    mc_PY[true_lep_idx],    mc_PZ[true_lep_idx],
                    mc_E[true_nu_idx],     mc_PX[true_nu_idx],     mc_PY[true_nu_idx],     mc_PZ[true_nu_idx],
                    mc_E[true_quark1_idx], mc_PX[true_quark1_idx], mc_PY[true_quark1_idx], mc_PZ[true_quark1_idx],
                    mc_E[true_quark2_idx], mc_PX[true_quark2_idx], mc_PY[true_quark2_idx], mc_PZ[true_quark2_idx],
                    }})
                    """)
        for name, omw in self._omega_wrappers.items():
            self.define_only_on(self._signal_categories, f"mc_sqme_{name}", omw, ["mc_ME_momenta", "mc_ME_flv"])
            # divide by recalculated nominal as all the ILD values are broken...
            if not name == "nominal":
                # TODO: get sqme from whizard
                self.define_only_on(self._signal_categories, f"weight_{name}", f"mc_sqme_{name} / mc_sqme_nominal")
