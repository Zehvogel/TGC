import ROOT
import numpy as np
import scipy as sc
from array import array


def make_THnSparse_like(h):
    axes = list(h.GetListOfAxes())
    return ROOT.THnSparseD("", "", axes)


def make_observed_histogram(histos: dict, meta: dict, run_config: dict[str, float]):
    """Merges histograms for 4 helicities into one according to a run_config containing lumi and beam pols"""
    # get axes from first histogram
    example_hist = list(histos.values())[0]
    # make empty sparse hist
    res = make_THnSparse_like(example_hist)
    for name, h in histos.items():
        m = meta[name]
        lumi_weight =  run_config["lumi"] / m["lumi"].GetVal()
        pol_weight = 0.25 * (1.0 + run_config["e_pol"] * m["e_pol"].GetVal()) * (1.0 + run_config["p_pol"] * m["p_pol"].GetVal())
        weight = lumi_weight * pol_weight
        # OMG this even adds bins with content 0
        # https://github.com/root-project/root/issues/19366
        res.Add(h, weight)
    return res


def make_THn_like(h):
    axes = list(h.GetListOfAxes())
    nbins = array("i", [a.GetNbins() for a in axes])
    xmin = array("d", [a.GetXmin() for a in axes])
    xmax = array("d", [a.GetXmax() for a in axes])
    return ROOT.THnD("", "",len(axes), nbins, xmin, xmax)


def make_observed_histogram_fast(histos: dict, meta: dict, run_config: dict[str, float]):
    """Merges histograms for 4 helicities into one according to a run_config containing lumi and beam pols"""
    # get axes from first histogram
    example_hist = list(histos.values())[0]
    # make empty sparse hist
    res = make_THn_like(example_hist)
    for name, h in histos.items():
        m = meta[name]
        lumi_weight =  run_config["lumi"] / m["lumi"].GetVal()
        pol_weight = 0.25 * (1.0 + run_config["e_pol"] * m["e_pol"].GetVal()) * (1.0 + run_config["p_pol"] * m["p_pol"].GetVal())
        weight = lumi_weight * pol_weight
        res.Add(h, weight)
    return res


def hist_to_np(h):
    axes = h.GetListOfAxes()
    n = len(axes)
    print(axes, n)
    nbins = [a.GetNbins() + 2 for a in axes]
    # print(nbins)
    res = np.zeros(nbins)
    # print(h.GetNbins())
    for i in range(h.GetNbins()):
        bc = h.GetBinContent(i, 0)
        res.flat[i] = bc

    return res


def make_asimov_np(h_np):
    return sc.stats.poisson.rvs(h_np, size=h_np.shape)


def make_asimov(h, seed: int):
    rnd = ROOT.TRandomMT64(seed)
    res = make_THnSparse_like(h)
    coord = array("i", [0] * h.GetNdimensions())
    for i in range(h.GetNbins()):
        val = h.GetBinContent(i, coord)
        new_val = rnd.Poisson(val)
        if new_val != 0:
            res.SetBinContent(coord, new_val)
    res.Scale(h.Integral(True) / res.Integral(True))
    return res


def make_asimov_fast(h, seed: int):
    rnd = ROOT.TRandomMT64(seed)
    res = make_THn_like(h)
    for i in range(h.GetNbins()):
        val = h.GetBinContent(i, 0)
        new_val = rnd.Poisson(val)
        if new_val != 0:
            res.SetBinContent(i, new_val)
    res.Scale(h.Integral(True) / res.Integral(True))
    return res


def make_1D_projections(hn):
    n = hn.GetNdimensions()
    res = []
    for i in range(n):
        h = hn.Projection(i)
        h.SetDirectory(0)
        res.append(h)
    return res


def make_1d_projections_np(hn_np):
    n = hn_np.ndim
    res = []
    for i in range(n):
        proj = np.sum(hn_np, axis=tuple([j for j in range(n) if j != i]))
        # cut off flow bins
        res.append(proj[1:-1])
    return res


def make_observables_np(h_np, bin_centers):
    res = []
    for h in h_np:
        res.append(np.dot(h, bin_centers))
    return res

def make_new_asimov_observables_np(h, bin_centers):
    h_asimov = make_asimov_np(h)
    h_1Ds = make_1d_projections_np(h_asimov)
    return make_observables_np(h_1Ds, bin_centers)


def make_observable(h):
    res = 0.
    for i in range(h.GetNbinsX()):
        bin_content = h.GetBinContent(i+1)
        bin_center = h.GetBinCenter(i+1)
        # print(f"i: {i} bin_center: {bin_center}, bin_content: {bin_content}")
        res += bin_center * bin_content
    return res

def make_observables(histos):
    res = []
    for h in histos:
        res.append(make_observable(h))
    return res

def make_new_asimov_observables(h, seed: int):
    h_asimov = make_asimov(h, seed)
    h_1Ds = make_1D_projections(h_asimov)
    return make_observables(h_1Ds)

