# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 16:54:09 2025

@author: bbuky

This file contains the minimum amount of code required to analyse LBT signals using SPRINT.
"""
import matplotlib.pyplot as plt
from LBT_analyser import LBT_analyser
from build_LBT import BB_file_picker
#%% read parameter file
from parameterFile_SOUL_I_Band import initializeParameterFile
param = initializeParameterFile()
#%% Create your LBT

# This determines which telescope you're using and what binning you want
binning=1
BB_file_picker(param, 'Left', binning)

# This initializes the class and creates the desired LBT model
LBT = LBT_analyser(param,binning,make_plots=True,atm=False)

#%% Initialize SPRINT

LBT.init_SPRINT(recompute_sensitivity=True)

#%% Use SPRINT to estimate mis-registrations for raw signals from LBT

# Define the files to use
trck = '20201001_075153'

loc = '../demodulated_slopes_orig/bin_'+str(LBT.binning)+'/'+str(trck)

slopes = loc + '/demodulated_slopes.fits'
phi = loc + '/phi.fits'
info = loc + '/data_info.fits'

# De-modulate the LBT data
LBT.get_on_sky_modulated_signal(slopes=slopes, phi=phi, info=info)

# Plot the signal being sent to SPRINT
plt.figure()
plt.imshow(LBT.slopes_2D)
plt.title('Signal')
plt.show()

# Run SPRINT
LBT.run_SPRINT(n_iteration=6,gain_estimation=1)

print('Flux Method shift estimates = ' + str(round(LBT.data_info.sx,4)) + ',' + str(round(LBT.data_info.sy,4)))
print('SPRINT estimate (X shift, Y shift, rotation) = ' + str(LBT.misreg_est))
