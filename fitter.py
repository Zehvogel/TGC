import ROOT
from array import array
from math import sqrt


class FitHandler:
    _signal_histograms_nd: dict
    _template_parametrisations: dict
    _background_histograms_nd: dict
    _oo_matrices: dict
    n_couplings: float
    n_bins: float
    obs_names: list[str]
    _signal_meta: dict
    parameters: dict
    process_e_pol = [-1., -1., 1., 1.]
    process_p_pol = [-1., 1., -1., 1.]

    def __init__(self, file_path: str, config: dict):
        ROOT.gInterpreter.Declare("#include \"fit.h\"")
        self._signal_histograms_nd = {}
        self._background_histograms_nd = {}
        self._oo_matrices = {}
        self._template_parametrisations = {}
        self.n_couplings = len(config["parameters"])
        self.n_bins = config["n_bins"]
        self.obs_names = config["obs_names"]
        self._signal_meta = {}
        self.parameters = config["parameters"]
        with ROOT.TFile(file_path) as input_file:
            # take care of signals first
            signal_dir = input_file[config["signal_cat"]]
            for process_name in config["signal_processes"]:
                p_dir = signal_dir[process_name]
                obs = "obs_ND"
                signal_histo = p_dir[obs]
                meta_dir = p_dir["meta"]
                obs_meta = {}
                for key in meta_dir.GetListOfKeys():
                    key_name = key.GetName()
                    obs_meta[key_name] = meta_dir[key_name]
                # normalise histograms to 1 fb_inv
                signal_histo.Scale(1 / obs_meta["lumi"].GetVal())
                self._signal_histograms_nd[process_name] = signal_histo
                self._signal_meta[process_name] = obs_meta
                # get OO matrix
                mat_dir = input_file["oo_matrix"]
                mat = mat_dir[process_name]
                self._oo_matrices[process_name] = mat / obs_meta["lumi"].GetVal()
            # get template parameters
            # unfortunately here we have the loops the other way around
            template_dir = input_file["template_parametrisations"]
            for obs in config["obs_names"]:
                obs_dir = template_dir[obs]
                p_pars = {}
                for process_name in config["signal_processes"]:
                    p_dir = obs_dir[process_name]
                    pars = {}
                    for par in config["parameters"].keys():
                        par_hist = p_dir[par]
                        # very consistent ownership model of root requires us to do this...
                        par_hist.SetDirectory(ROOT.nullptr)
                        pars[par] = par_hist
                    p_pars[process_name] = pars
                self._template_parametrisations[obs] = p_pars
            for bkg in config["backgrounds"]:
                bkg_dir = input_file[bkg]
                bkg_hist = [self.make_THn_like(list(self._signal_histograms_nd.values())[0]) for i in range(4)]
                for process in bkg_dir.GetListOfKeys():
                    p_dir = bkg_dir[process.GetName()]
                    # should contain meta, O_bla... and obs_ND
                    meta = p_dir["meta"]
                    h = p_dir["obs_ND"]
                    # normalise histograms to 1 fb_inv
                    h.Scale(1. / meta["lumi"].GetVal())
                    # figure out which pol hist to add this to from e_pol and_p_pol
                    e_pol = meta["e_pol"].GetVal()
                    p_pol = meta["p_pol"].GetVal()
                    bkg_hist[self.pol_idx(e_pol, p_pol)].Add(h)
                self._background_histograms_nd[bkg] = bkg_hist


    @staticmethod
    def make_THn_like(h):
        axes = list(h.GetListOfAxes())
        nbins = array("i", [a.GetNbins() for a in axes])
        xmin = array("d", [a.GetXmin() for a in axes])
        xmax = array("d", [a.GetXmax() for a in axes])
        return ROOT.THnD("", "",len(axes), nbins, xmin, xmax)


    @staticmethod
    def pol_idx(e_pol: float, p_pol: float):
        return int(((e_pol + 1) * 2 + (p_pol + 1))/2)


    # def make_expected_fun(self, obs_name: str, run_config: dict[str, float]):
    def make_expected_fun(self, obs_name: str):
        signal_hists = self.get_signal_hist_1d(obs_name)
        template_pars = self.get_template_params(obs_name)
        # lumi = run_config["lumi"]
        # e_pol = run_config["e_pol"]
        # p_pol = run_config["p_pol"]
        n_parameters = self.n_couplings + 3 + 1
        # fun = ROOT.expect_fun[f"{self.n_bins}, {self.n_couplings}"](signal_hists, template_pars, lumi, e_pol, p_pol)
        fun = ROOT.expect_fun[f"{self.n_bins}, {self.n_couplings}"](signal_hists, template_pars)
        functor = ROOT.Math.Functor(fun, n_parameters)
        return fun, functor


    def get_signal_hist_1d(self, obs_name: str):
        obs_idx = self.obs_names.index(obs_name)
        projections = []
        for h in self._signal_histograms_nd.values():
            p = h.Projection(obs_idx)
            p.SetDirectory(0)
            projections.append(p)
        return projections


    def get_template_params(self, obs_name: str):
        tmp_pars = self._template_parametrisations[obs_name]
        # need to convert dict of dicts to list of lists
        return [list(p.values()) for p in tmp_pars.values()]


    def build_model(self, run_configs: list(dict[str, float])):

        # TODO: separate run components and build roo simultaneous by hand

        # build all the RooFit stuff
        # first build initial obs and the covariance matrix according to the run config
        obs_initial = self.get_initial_observables(run_config)
        cov_matrix = self.get_obs_cov_matrix(run_config)

        # define observables
        obs_pars = []
        for i, o in enumerate(obs_initial):
            name = self.obs_names[i]
            sigma = sqrt(cov_matrix[i][i])
            min_o = o - 3 * sigma
            max_o = o + 3 * sigma
            ob = ROOT.RooRealVar(name, name, min_o, max_o)
            obs_pars.append(ob)

        # self.obs_pars = obs_pars

        # define nuisance parameters for run parameters
        lumi_par = ROOT.RooRealVar("lumi", "lumi", run_config["lumi"], 0.9 * run_config["lumi"], 1.1 * run_config["lumi"])
        lumi_par.setConstant()
        e_pol_par = ROOT.RooRealVar("e_pol", "e_pol", run_config["e_pol"], -1., 1.)
        e_pol_par.setConstant()
        p_pol_par = ROOT.RooRealVar("p_pol", "p_pol", run_config["p_pol"], -1., 1.)
        p_pol_par.setConstant()

        # define coupling parameters
        coupling_pars = []
        for name, value in self.parameters.items():
            # that should be tight enough...
            cpl = ROOT.RooRealVar(name, name, value, -0.5, 0.5)
            coupling_pars.append(cpl)
        # self.coupling_pars = coupling_pars

        # define nuisance parameters
        # TODO rename because these are normfactors?
        nuisance_pars = []
        # signal only so far
        nu_par = ROOT.RooRealVar("mu_signal", "mu_signal", 1., 0.9, 1.1)
        nu_par.setConstant()
        nuisance_pars.append(nu_par)
        all_pars = [lumi_par, e_pol_par, p_pol_par] + coupling_pars + nuisance_pars
        # self.obs_and_pars = obs_pars + all_pars

        self.expected_funs = []
        bound_expected_funs = []
        for i, name in enumerate(self.obs_names):
            fun, func = self.make_expected_fun(name, run_config)
            self.expected_funs.append(fun)
            self.expected_funs.append(func)
            bound_fun = ROOT.RooFit.bindFunction(f"expectation_{name}", func, all_pars)
            # self.expected_funs.append(bound_fun)
            bound_expected_funs.append(bound_fun)

        print("hello")
        gauss = ROOT.RooMultiVarGaussian("multi_gauss", "multi_gauss", obs_pars, bound_expected_funs, cov_matrix)
        # self.gauss = gauss

        # global obs for pol
        # g_obs = []
        # e_pol_glob = ROOT.RooRealVar("e_pol_glob", "e_pol_glob", run_config["e_pol"], -1., 1.)
        # e_pol_glob.setConstant()
        # g_obs.append(e_pol_glob)
        # e_pol_const = ROOT.RooGaussian("e_pol_const", "e_pol_const", e_pol_glob, e_pol_par, 0.1)

        # constraints = []
        # constraints.append(e_pol_const)

        # self.g_obs = g_obs
        # self.constraints = constraints

        # model = ROOT.RooProdPdf("model", "model", [gauss] + constraints)
        model = gauss

        ws = ROOT.RooWorkspace("ws")
        ws.Import(model)
        # TODO: model config

        return ws


    def get_obs_cov_matrix(self, run_config):
        n_obs = len(self.obs_names)
        sum_mat = ROOT.Math.SVector[f"double, {int(n_obs * (n_obs + 1) / 2)}"]()
        for i, mat in enumerate(self._oo_matrices.values()):
            e_pol = self.process_e_pol[i]
            p_pol = self.process_p_pol[i]
            pol_weight = 0.25 * (1.0 + run_config["e_pol"] * e_pol) * (1.0 + run_config["p_pol"] * p_pol)
            # just need to multiply by lumi as they are already normalised
            weight = run_config["lumi"] * pol_weight
            mat *= weight
            sum_mat += mat
        # first build smatrix, then convert it to tmatrix...
        smat = ROOT.Math.SMatrix[f"double, {n_obs}, {n_obs}, ROOT::Math::MatRepSym<double, {n_obs}"](sum_mat, False)
        tmat = ROOT.TMatrixDSym(n_obs)
        for i in range(n_obs):
            for j in range(n_obs):
                tmat[i][j] = smat[i][j]
        return tmat


    def get_initial_observables(self, run_config: dict[str, float]):
        # TODO background
        obs_hist = FitHandler.make_observed_histogram_fast(self._signal_histograms_nd, self._signal_meta, run_config)
        h_1d = FitHandler.make_1D_projections(obs_hist)
        return FitHandler.make_observables(h_1d)


    @staticmethod
    def make_observed_histogram_fast(histos: dict, meta: dict, run_config: dict[str, float]):
        """Merges histograms for 4 helicities into one according to a run_config containing lumi and beam pols"""
        # get axes from first histogram
        example_hist = list(histos.values())[0]
        # make empty hist
        res = FitHandler.make_THn_like(example_hist)
        for name, h in histos.items():
            m = meta[name]
            pol_weight = 0.25 * (1.0 + run_config["e_pol"] * m["e_pol"].GetVal()) * (1.0 + run_config["p_pol"] * m["p_pol"].GetVal())
            # just need to multiply by lumi as they are already normalised
            weight = run_config["lumi"] * pol_weight
            res.Add(h, weight)
        return res


    @staticmethod
    def make_observable(h):
        res = 0.
        for i in range(h.GetNbinsX()):
            bin_content = h.GetBinContent(i+1)
            bin_center = h.GetBinCenter(i+1)
            # print(f"i: {i} bin_center: {bin_center}, bin_content: {bin_content}")
            res += bin_center * bin_content
        return res


    @staticmethod
    def make_observables(histos):
        res = []
        for h in histos:
            res.append(FitHandler.make_observable(h))
        return res


    @staticmethod
    def make_1D_projections(hn):
        n = hn.GetNdimensions()
        res = []
        for i in range(n):
            h = hn.Projection(i)
            h.SetDirectory(0)
            res.append(h)
        return res