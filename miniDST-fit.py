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


model = fit_handler.build_model(run)


# In[ ]:


# model.fitTo(ds)
# model.getVal()
model.Print("t")
