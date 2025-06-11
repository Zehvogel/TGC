from analysis_framework.Analysis import Analysis
import ROOT
import numpy as np
from OO.whizard.model_parser import ModelParser
import subprocess
from alt_setup_creator import AltSetupHandler


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
    _template_graphs = {}
    _template_funs = {}

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


    def define_truth_objects(self, categories: list[str]):
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


        def first_two_unstable_between(abs_pdg_min, abs_pdg_max):
            return f"""
            auto& genStat = MCParticlesSkimmed.generatorStatus;
            auto& pdg = MCParticlesSkimmed.PDG;
            auto mask = genStat == 2 && abs(pdg) < {abs_pdg_max} && abs(pdg) > {abs_pdg_min};
            // abuse ArgMax to get the first set position
            auto idx = ArgMax(mask);
            // Drop the first index and use ArgMax again
            auto idx2 = ArgMax(Drop(mask, {{idx}}));
            // increment by one to compensate for removal of the first
            idx2++;
            return ROOT::VecOps::RVec({{idx, idx2}});
            """


        self._define(("true_lep_idcs", first_two_unstable_between(10, 13)), categories)
        self._define(("true_lep_idx", "true_lep_idcs[0]"), categories)
        self._define(("true_nu_idx", "true_lep_idcs[1]"), categories)
        self._define(("true_quarks_idcs", first_two_unstable_below(6)), categories)
        self._define(("true_quark1_idx", "true_quarks_idcs[0]"), categories)
        self._define(("true_quark2_idx", "true_quarks_idcs[1]"), categories)
        self._define(("true_lep_lvec", make_lvec_M("MCParticlesSkimmed", "true_lep_idx")), categories)
        self._define(("true_nu_lvec", make_lvec_M("MCParticlesSkimmed", "true_nu_idx")), categories)
        self._define(("true_quark1_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark1_idx")), categories)
        self._define(("true_quark2_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark2_idx")), categories)
        self._define(("true_leptonic_W_lvec", "true_lep_lvec + true_nu_lvec"), categories)
        self._define(("true_hadronic_W_lvec", "true_quark1_lvec + true_quark2_lvec"), categories)
        self._define(("true_iso_lep_charge", "MCParticlesSkimmed.PDG[true_lep_idx] > 0. ? -1. : 1."), categories)

        self.truth_defined = True
        self.truth_categories = categories


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
            # if not name == "nominal":
            # divide by recalculated nominal as all the ILD values are broken...
            self.define_only_on(self._signal_categories, f"weight_{name}", f"mc_sqme_{name} / mc_sqme_nominal")


    def define_optimal_observables(self, alt_handler: AltSetupHandler, only=None):
        for name, var in alt_handler.get_named_variations_it(only):
            self.Define(f"O_{name}", f"{1/var} * (reco_sqme_12_nominal + reco_sqme_21_nominal - reco_sqme_12_{name} - reco_sqme_21_{name}) / (reco_sqme_12_nominal + reco_sqme_21_nominal)")
            self._define((f"mc_O_{name}", f"{1/var} * (mc_sqme_nominal - mc_sqme_{name}) / mc_sqme_nominal"), self._signal_categories)


    def plot_template_bins(self, observable_name: str, alt_handler: AltSetupHandler, only=None):
        template_name = f"tmpl_{observable_name}"
        varied_histos = self._varied_histograms[template_name]
        canvases = {}
        for process_name, v_histos in varied_histos.items():
            nominal_histo = v_histos["nominal"]
            if only:
                raise NotImplementedError
            par_graphs = {}
            par_canvases = {}
            par_funs = {}
            for par_name, vars in alt_handler.get_variations_ext().items():
                x_vals = sorted(vars)
                r_histos = []
                for v in x_vals:
                    var_name = alt_handler.make_name(par_name, v)
                    r_histos.append(v_histos[f"weight:{var_name}"] / nominal_histo)

                bin_graphs = []
                bin_canvases = []
                bin_funs = []
                for i in range(1, nominal_histo.GetNcells()-1):
                    graph = ROOT.TGraph()
                    for j, v in enumerate(x_vals):
                        if j > 0 and v == -x_vals[j-1]:
                            graph.AddPoint(0., 1.)
                        graph.AddPoint(v, r_histos[j].GetBinContent(i))
                    fun = ROOT.TF1("", "1 + [0]*x + [1]*x*x", x_vals[0], x_vals[-1])
                    fun.SetParameter(0, 1.)
                    fun.SetParameter(1, 0.1)
                    graph.Fit(fun)
                    bin_funs.append(fun)
                    bin_graphs.append(graph)
                    c = ROOT.TCanvas()
                    graph.Draw("apl")
                    # fun.Draw("same")
                    graph.SetTitle(f"{process_name} bin {i};{par_name}")
                    c.Draw()
                    bin_canvases.append(c)
                    # print("hello")
                par_graphs[par_name] = bin_graphs
                par_canvases[par_name] = bin_canvases
                par_funs[par_name] = bin_funs
            self._template_graphs[process_name] = par_graphs
            canvases[process_name] = par_canvases
            self._template_funs[process_name] = par_funs
        self._canvases["template_bins"] = canvases

