from analysis_framework.Analysis import Analysis
import ROOT
import numpy as np
from OO.whizard.model_parser import ModelParser
import subprocess
from alt_setup_creator_new import AltSetupHandler
from itertools import combinations_with_replacement, product
from math import sqrt
from array import array
from time import sleep

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
    _template_pars = {}
    _oo_matrix = {}
    _oo_matrix_var = {}

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


    def define_optimal_observables(self, names: str):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self.Define(f"O_{name}", f"{1/var} * (reco_sqme_12_nominal + reco_sqme_21_nominal - reco_sqme_12_{name} - reco_sqme_21_{name}) / (reco_sqme_12_nominal + reco_sqme_21_nominal)")
            self._define((f"mc_O_{name}", f"{1/var} * (mc_sqme_nominal - mc_sqme_{name}) / mc_sqme_nominal"), self._signal_categories)


    def book_oo_matrix(self, observables: list[str]):
        # TODO: all or only signal???
        # FIXME: well this would be the clean way but all the symmetric matrix implementations in ROOT suck big time
        elements= ["*".join(c) for c in combinations_with_replacement(observables, 2)]
        # elements= ["*".join(c) for c in product(observables, repeat=2)]
        self._define(("oo_matrix", f"ROOT::RVecD{{{','.join(elements)}}}"), self._signal_categories)
        n = len(elements)
        res = ROOT.RVecD(n, 0.)
        # f = ROOT.VecOps.__add__
        # f = lambda x, y: x + y
        # ROOT.gInterpreter.Declare("auto foo = [] (const ROOT::RVecD& a, const ROOT::RVecD& b) {return a+b;};")
        # ROOT.gInterpreter.Declare("""
        # ROOT::RVecD rvec_reduce_d(ROOT::RVecD a, ROOT::RVecD b) {
        #     return a + b;
        # }
        # """)
        ROOT.gInterpreter.Declare("""
        struct rvec_reducer_d {
            ROOT::RVecD operator()(ROOT::RVecD a, ROOT::RVecD b) {
            return a + b;
            }
        };
        """)
        f = ROOT.rvec_reducer_d()
        self._oo_matrix = self.book_some_method("Reduce", (f, "oo_matrix", res), categories=self._signal_categories)
        # also make linear parametrisations
        for o in observables:
            var_name = o.removeprefix("O_")
            self._define((f"oo_matrix_{var_name}", f"oo_matrix * weight_{var_name}"), self._signal_categories)
            self._oo_matrix_var[var_name] = self.book_some_method("Reduce", (f, f"oo_matrix_{var_name}", res), categories=self._signal_categories)


    def calculate_template_parametrisation(self, observable_name: str, alt_handler: AltSetupHandler):
        template_name = f"tmpl_{observable_name}"
        varied_histos = self._varied_histograms[template_name]
        template_graphs = {}
        template_funs = {}
        template_pars = {}
        for process_name, v_histos in varied_histos.items():
            nominal_histo = v_histos["nominal"]
            par_graphs = {}
            par_funs = {}
            par_pars = {}
            for par_name, vars in alt_handler.get_variations_ext().items():
                par_histo = nominal_histo.Clone()
                x_vals = sorted(vars)
                r_histos = []
                for v in x_vals:
                    var_name = alt_handler.make_name(par_name, v)
                    r_histos.append(v_histos[f"weight:{var_name}"] / nominal_histo)

                bin_graphs = []
                bin_funs = []
                # ignore over/underflow bins
                n_bins = nominal_histo.GetNcells() - 2
                for i in range(1, n_bins + 1):
                    graph = ROOT.TGraph()
                    for j, v in enumerate(x_vals):
                        if j > 0 and v == -x_vals[j-1]:
                            graph.AddPoint(0., 1.)
                        graph.AddPoint(v, r_histos[j].GetBinContent(i))
                    fun = ROOT.TF1("", "1 + [0]*x + [1]*x*x", x_vals[0], x_vals[-1])
                    fun.SetParameter(0, 1.)
                    fun.SetParameter(1, 0.1)
                    graph.Fit(fun, "Q")
                    par1 = fun.GetParameter(0)
                    par2 = fun.GetParameter(1)
                    par_histo.SetBinContent(i, par1)
                    # print(f"parameters after fit: ({par1}, {par2})")
                    bin_funs.append(fun)
                    bin_graphs.append(graph)
                par_graphs[par_name] = bin_graphs
                par_funs[par_name] = bin_funs
                par_pars[par_name] = par_histo
            template_graphs[process_name] = par_graphs
            template_funs[process_name] = par_funs
            template_pars[process_name] = par_pars
        self._template_graphs[observable_name] = template_graphs
        self._template_funs[observable_name] = template_funs
        self._template_pars[observable_name] = template_pars



    def calculate_template_parametrisation_full(self, observable_name: str, alt_handler: AltSetupHandler):
        template_name = f"tmpl_{observable_name}"
        varied_histos = self._varied_histograms[template_name]
        pars = alt_handler.get_pars()
        n_pars = len(pars)

        # urgh no point in doing this as ROOT sorts the parameters internally anyway back into the same order as with the other definition.......
        terms = []
        for k, (par1, par2) in enumerate(combinations_with_replacement(["1"] + alt_handler.get_pars(), r=2)):
            # print(k, par1, par2)
            if k == 0:
                continue
            if par1 == "1":
                i = pars.index(par2)
                # terms.append(f"x[{i}] * [par_{i}]")
                terms.append(f"x[{i}] * [p_{k-1}]")
            else:
                i = pars.index(par1)
                j = pars.index(par2)
                # terms.append(f"x[{i}] * x[{j}] * [par_{i}_{j}]")
                terms.append(f"x[{i}] * x[{j}] * [p_{k-1}]")
        parametrisation = " + ".join(terms)

        # parametrisation = ""
        # for i in range(n_pars):
        #     if i > 0:
        #         parametrisation += " + "
        #     parametrisation += f"x[{i}] * [par_{i}]"
        #     for j in range(n_pars):
        #         if j < i:
        #             continue
        #         parametrisation += f" + x[{i}] * x[{j}] * [par_{i}_{j}]"
        print(parametrisation)
        # un-documented root feature, set min > max and more than 1D is fine...
        # otherwise there is an error and the zombie bit gets set
        f = ROOT.TF1("f", parametrisation, 1, 0)
        fit_fun = ROOT.Math.WrappedMultiTF1(f, f.GetNdim())
        fitter = ROOT.Fit.Fitter()
        fitter.SetFunction(fit_fun, False)
        template_pars = {}
        for process_name, v_histos in varied_histos.items():
            nominal_histo = v_histos["nominal"]
            fit_data = []
            n_bins = nominal_histo.GetNbinsX()
            n_variations = alt_handler.get_n_variations()
            for i in range(n_bins):
                fd = ROOT.Fit.BinData(n_variations +
                1, n_pars, ROOT.Fit.BinData.kNoError)
                fd.Add(array("d", [0.]*n_pars), 1.)
                fit_data.append(fd)
            for i, (var_name, vars) in enumerate(alt_handler.get_variations_nd()):
            # for each variation get the histogram
                h_var = v_histos[f"weight:{var_name}"] / nominal_histo
                for b in range(n_bins):
                    bc = h_var[b]
                    fit_data[b].Add(array("d", vars), bc)
            fit_results = []
            for b in range(n_bins):
                fit_worked = fitter.Fit(fit_data[b])
                if not fit_worked:
                    print("ERROR fit returned false!!!")
                    # print(f"fit res: ")
                fitter.Result().Print(ROOT.std.cout)
                fit_result = fitter.Result()
                fit_results.append(list(fit_result.Parameters()))
            process_par_histos = {}
            print(fit_results)
            for i, (par1, par2) in enumerate(combinations_with_replacement(["1"] + alt_handler.get_pars(), r=2)):
                    # print(i, par1, par2)
                    if i == 0:
                        # nothing to do for no variation
                        continue
                    var_idx = i - 1
                    par_histo = nominal_histo.Clone()
                    for b in range(n_bins):
                        bc = fit_results[b][var_idx]
                        par_histo.SetBinContent(b+1, bc)
                    name = par2 if par1 == "1" else f"{par1}_{par2}"
                    par_histo.SetTitle(name)
                    process_par_histos[name] = par_histo
            template_pars[process_name] = process_par_histos
                    # fit
                    # get parameters from fit result
                    # put the parameters directly into some kind of matrix?
                    # or histograms so that we can view them in the file?
        self._template_pars[observable_name] = template_pars


    def plot_template_bins(self, observable_name: str, plot_path: str|None = None):
        canvases = {}
        template_graphs = self._template_graphs[observable_name]
        for process_name, par_graphs in template_graphs.items():
            par_canvases = {}
            for par_name, bin_graphs in par_graphs.items():
                n_bins = len(bin_graphs)
                # x_width = min((n_bins+1)*ROOT.gStyle.GetCanvasDefW() // 2, 8000)
                x_width = (n_bins+1)*ROOT.gStyle.GetCanvasDefW() // 2
                bin_canvas = ROOT.TCanvas("", "", x_width, ROOT.gStyle.GetCanvasDefH())
                bin_canvas.Divide(n_bins+1, 1, 0, 0)
                for i in range(1, n_bins + 1):
                    graph = bin_graphs[i-1]
                    bin_canvas.cd(i + 1)
                    graph.Draw("apl")
                    graph.SetTitle(f"{process_name} bin {i};{par_name}")
                # TODO: scale all graphs to the same y-axis and plot it on the first canvas...
                bin_canvas.Draw()
                if plot_path:
                    bin_canvas.SaveAs(f"{plot_path}/{observable_name}_{process_name}_{par_name}.pdf")
                par_canvases[par_name] = bin_canvas
            canvases[process_name] = par_canvases
        self._canvases[f"{observable_name}_template_bins"] = canvases


    def write_fit_inputs(self, observable_names: list[str], output_path: str):
        self.store_raw_histograms(observable_names, output_path)
        # TODO: also write OO cov_matrix
        with ROOT.TFile(f"{output_path}/raw_histograms.root", "update") as output_file:
            dir = output_file.mkdir("template_parametrisations")
            for observable in observable_names:
                obs_dir = dir.mkdir(observable)
                for process_name, template_pars in self._template_pars[observable].items():
                    process_dir = obs_dir.mkdir(process_name)
                    process_dir.cd()
                    for name, h in template_pars.items():
                        h.Write(name)
            oo_mat_dir = output_file.mkdir("oo_matrix")
            oo_mat_dir.cd()
            for name, oo_mat in self._oo_matrix.items():
                # mat_par = ROOT.TParameter["ROOT::RVecD"](name, oo_mat.GetValue())
                # mat_par.Write()
                # smat = ROOT.Math.SMatrix["double", oo_mat.size()](oo_mat.begin(), oo_mat.end())
                # smat.Write(name)
                # mat = ROOT.TMatrixDSym(int(sqrt(oo_mat.size())), oo_mat.data())
                # mat.Write(name)
                mat = ROOT.Math.SVector[f"double, {oo_mat.size()}"](oo_mat.begin(), oo_mat.end())
                oo_mat_dir.WriteObject(mat, name)
            oo_mat_var_dir = output_file.mkdir("oo_matrix_var")
            for var_name, vars in self._oo_matrix_var.items():
                var_dir = oo_mat_var_dir.mkdir(var_name)
                var_dir.cd()
                for mat_name, oo_mat in vars.items():
                    mat = ROOT.Math.SVector[f"double, {oo_mat.size()}"](oo_mat.begin(), oo_mat.end())
                    var_dir.WriteObject(mat, mat_name)


