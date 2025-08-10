#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
from fitter import FitHandler


# In[2]:


# TODO: make it possible to fit multiple runs simultaneously
run = {
    "lumi": 5000,
    # "lumi": 25000,
    "e_pol": 0.,
    "p_pol": 0.,
}

# https://arxiv.org/pdf/1506.07830
# However, https://arxiv.org/pdf/2503.19983 cites the above but assigns 0.45 to both mixed pols
# LCVision scenario uses 3 ab_inv instead of 2
ilc_250_h20_lumi = 2000
ilc_250_h20 = [
    {
        "lumi": ilc_250_h20_lumi * 0.675,
        "e_pol": -0.8,
        "p_pol": 0.3,
    },
    {
        "lumi": ilc_250_h20_lumi * 0.225,
        "e_pol": 0.8,
        "p_pol": -0.3,
    },
    {
        "lumi": ilc_250_h20_lumi * 0.05,
        "e_pol": -0.8,
        "p_pol": -0.3,
    },
    {
        "lumi": ilc_250_h20_lumi * 0.05,
        "e_pol": 0.8,
        "p_pol": 0.3,
    },
]

# run = ilc_250_h20[0]


# In[3]:


conf = {
    "parameters": {
        "g1z": 0.0,
        "ka": 0.0,
        "la": 0.0,
    },
    "obs_names": [
        "O_g1z_pos_1em05",
        "O_ka_pos_1em05",
        "O_la_pos_1em05",
    ],
    "signal_cat": "4f_sw_sl_signal",
    "signal_processes": [
        "4f_sw_sl_eLpL_signal",
        "4f_sw_sl_eLpR_signal",
        "4f_sw_sl_eRpL_signal",
        "4f_sw_sl_eRpR_signal",
    ],
    "backgrounds": [
        "4f_sl_bkg",
        "4f_not_sl",
        "2f",
        "3f",
        "5f",
        "6f",
        "higgs",
    ],
    "n_bins": 65,
}
input_path = "data/histograms/full/raw_histograms.root"


# In[4]:


fit_handler = FitHandler(input_path, conf)


# In[ ]:


ws = fit_handler.build_model([run])
# coupling_pars = fit_handler.coupling_pars
coupling_pars = [ws.var(name) for name in conf["parameters"]]
obs_pars = [ws.var(name) for name in conf["obs_names"]]


# In[ ]:


ws.Print("t")


# In[ ]:


model = ws.pdf("multi_gauss")


# In[ ]:


ds = ROOT.RooStats.AsymptoticCalculator.GenerateAsimovData(model, obs_pars)
ds.Print("v")


# In[ ]:


fit_res = model.fitTo(ds, Save=True)


# In[ ]:


fit_res.covarianceMatrix().Print()


# In[ ]:


nll = model.createNLL(ds, EvalBackend="cpu")


# In[ ]:


nll_minimizer = ROOT.RooMinimizer(nll)


# In[ ]:


get_ipython().run_line_magic('time', '')
nll_minimizer.migrad()


# In[ ]:


get_ipython().run_line_magic('time', '')
# nll_minimizer.hesse()
# nll_minimizer.minos(coupling_pars)


# In[ ]:


# nll.Print("t")


# In[ ]:


pll0 = nll.createProfile({coupling_pars[0]})
pll1 = nll.createProfile({coupling_pars[1]})
pll2 = nll.createProfile({coupling_pars[2]})


# In[ ]:


frame0 = coupling_pars[0].frame(Range=(-0.004, 0.004))
frame1 = coupling_pars[1].frame(Range=(-0.004, 0.004))
frame2 = coupling_pars[2].frame(Range=(-0.004, 0.004))
nll.plotOn(frame0, ShiftToZero=True)
nll.plotOn(frame1, ShiftToZero=True)
nll.plotOn(frame2, ShiftToZero=True)


# In[ ]:


pll0.plotOn(frame0, LineColor="r")
pll1.plotOn(frame1, LineColor="r")
pll2.plotOn(frame2, LineColor="r")


# In[ ]:


c0 = ROOT.TCanvas()
frame0.SetMinimum(0)
frame0.SetMaximum(4)
frame0.Draw()
c0.Draw()
# c0.SaveAs("plots/fit/ll_pll.pdf(")
# c0.SaveAs("plots/fit/ll_pll_g1z.pdf")

c1 = ROOT.TCanvas()
frame1.SetMinimum(0)
frame1.SetMaximum(4)
frame1.Draw()
c1.Draw()
# c1.SaveAs("plots/fit/ll_pll.pdf")
# c1.SaveAs("plots/fit/ll_pll_ka.pdf")

c2 = ROOT.TCanvas()
frame2.SetMinimum(0)
frame2.SetMaximum(4)
frame2.Draw()
c2.Draw()
# c2.SaveAs("plots/fit/ll_pll.pdf)")
# c2.SaveAs("plots/fit/ll_pll_la.pdf")


# In[ ]:




