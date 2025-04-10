# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 16:54:09 2025

@author: bbuky

This file contains the minimum amount of code required to analyse LBT signals using SPRINT.
"""
import matplotlib.pyplot as plt
from LBT_analyser import LBT_analyser
from build_LBT import BB_file_picker
import numpy as np
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

#%% Use SPRINT to estimate mis-registrations of one set of raw signals from LBT

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

#%% Run SPRINT on a series of files (eg if a ramp of misregistrations have been applied)

# Define the files using the unique tracking number for each set
trck = ['20201001_075009','20201001_075108','20201001_075153','20201001_075233','20201001_075444']

# Create arrays for storing flux method and SPRINT estimates
Flux_ests = np.zeros((len(trck),2))
SPRINT_ests = np.zeros((len(trck),3))

for i in range(len(trck)):
    
    loc = '../demodulated_slopes_orig/bin_1/'+str(trck[i])
    
    slopes = loc + '/demodulated_slopes.fits'
    phi = loc + '/phi.fits'
    info = loc + '/data_info.fits'

    # De-modulate the LBT data
    LBT.get_on_sky_modulated_signal(slopes=slopes, phi=phi, info=info)

    # Plot the signal being sent to SPRINT
    plt.figure()
    plt.imshow(LBT.slopes_2D)
    plt.title('Signal for ' + str(trck[i]))
    plt.show()
    
    # Run SPRINT
    LBT.run_SPRINT(n_iteration=6,gain_estimation=1)
    print('SPRINT MISREGS = ' + str(LBT.misreg_est))
    
    Flux_ests[i] = [LBT.data_info.sx,LBT.data_info.sy]
    SPRINT_ests[i] = LBT.misreg_est
    
fig,axs = plt.subplots(1,2,figsize=(10,5))

axs[0].plot(SPRINT_ests[:,0],'.--',label='SPRINT')
axs[0].plot(Flux_ests[:,0],'.--',label='LBT')
axs[0].set(xlabel='Mis-registration case',ylabel='ShiftX [m]')
axs[0].legend()

axs[1].plot(SPRINT_ests[:,1],'.--',label='SPRINT')
axs[1].plot(Flux_ests[:,1],'.--',label='LBT')
axs[1].set(xlabel='Mis-registration case',ylabel='ShiftY [m]')
axs[1].legend()

plt.suptitle('SPRINT vs Flux Method on LBT Data')
fig.tight_layout()
plt.show()
    
