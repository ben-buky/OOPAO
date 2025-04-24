# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 14:30:43 2025

@author: tsl29789
"""

import matplotlib.pyplot as plt
from LBT_analyser import LBT_analyser
from build_LBT import BB_file_picker
import numpy as np
from OOPAO.MisRegistration import MisRegistration
from OOPAO.mis_registration_identification_algorithm.applyMisRegistration import applyMisRegistration
#%% read parameter file
from parameterFile_SOUL_I_Band import initializeParameterFile
param = initializeParameterFile()

#%% Create your LBT

# Set your binning, whether you're using new or old data, and which telescope you want 
# Then use BB_file_picker to generate the correct file paths parameters
param['new_IF'] = True
param['which_tel'] = 'Left'
binning=2
BB_file_picker(param, binning)

m_ref_og = MisRegistration()
m_ref_og.shiftX            = 0.137
m_ref_og.shiftY            = 0.023
m_ref_og.rotationAngle     = 209.47
m_ref_og.radialScaling     = 0.025  
m_ref_og.tangentialScaling = 0.019 

m_ref = MisRegistration()
m_ref.shiftX            = 0.137
m_ref.shiftY            = 0.023
m_ref.rotationAngle     = 209.47
m_ref.radialScaling     = 0.025  
m_ref.tangentialScaling = 0.019 

# This initializes the class and creates the desired LBT model
LBT_new = LBT_analyser(param,binning,make_plots=True,atm=False,misReg=m_ref)

LBT_new.m_ref.print_()

#%% Initialize SPRINT for reference if doing in full

LBT_new.init_SPRINT_ref(recompute_sensitivity=True)

#%% Run reference calculation

LBT_new.run_ref(n_iteration=6, n_update_zero_point=2)

#%% Initialize Sprint for several modes to do cruder but faster reference calculation
modes = [1,2,3,4,5,6,7,8,9,10,30,40,60,80,110,150,193,218]
LBT_new.init_SPRINT(mode=modes, n_mis_reg=3, recompute_sensitivity=True)

#%% 
from OOPAO.tools.displayTools import display_wfs_signals
# Plot the signal being sent to SPRINT in 2D
display_wfs_signals(LBT_new.wfs, signals = LBT_new.calib_CL.D[:,modes]),plt.colorbar()
display_wfs_signals(LBT_new.wfs, signals = LBT_new.sprint.calib_0.D),plt.colorbar()

#%% Run SPRINT on the reference signal with zeropoint updates to try and find reference

LBT_new.on_sky_slopes = LBT_new.calib_CL.D[:,modes]

# Run SPRINT - SPRINT will run on whatever is saved as LBT.on_sky_slopes
LBT_new.run_SPRINT(n_iteration=10, n_update_zero_point=0, gain_estimation=1)