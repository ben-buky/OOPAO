# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 16:54:09 2025

@author: tsl29789

This tutorial aims to explain in as much detail as possible how to use LBT_analyser to run SPRINT on some LBT data.
"""
import matplotlib.pyplot as plt
from LBT_analyser import LBT_analyser
from build_LBT import BB_file_picker
#%% read parameter file
from parameterFile_SOUL_I_Band import initializeParameterFile
param = initializeParameterFile()
#%% Create your LBT

# Set your binning and which telescope you want to generate the correct file paths parameters
binning=1
BB_file_picker(param, 'Left', binning)

""" If you have a different file structure to Ben, you will need to adapt BB_file_picker or define the extra parameters you need manually.
The required extra file path parameters are:
    
    param['filename_pup'] - pupil mask
    param['int_mat'] - interaction matrix
    param['slopex'] - also used for generating interaction matrix, specific to the binning used
    param['filename_if'] - influence functions
    param['filename_mir_modes'] - ASM eigenmodes
    param['filename_m2c'] - Mode to Command matrix
    param[filename_coord] - Coordinates of ASM actuators
    
"""

# This initializes the class and creates the desired LBT model
LBT = LBT_analyser(param=param,     # define where to find the model parameters
                   binning=binning, # set your binning, must match what you use for BB_file_picker
                   misReg=None,     # you can provide a unique starting misregistration, None will use the standard zero-point reference for your binning
                   psim=False,      # you have the option to generate a synthetic IM rather than using the real one, the default is False
                   make_plots=True, # you can turn off the plots produced in the setup
                   n_modes=500,     # specify your number of modes, this will default to the maximum number allowed for your binning
                   atm=True)        # you can prevent an atmosphere being generated if you're only running SPRINT and want to speed up the code

# The above initialization stores your parameters as LBT.param and binning as LBT.binning. 
# You can also now access objects of the model and the reference mis-registration through your LBT, some examples are below:

# Display the telescope OPD 
LBT.tel + LBT.atm # will only work if you've generated an atmosphere
LBT.ngs*LBT.tel
plt.figure()
plt.imshow(LBT.tel.OPD)
plt.title('LBT OPD with atm')
plt.show()

# Display the WFS camera

plt.figure()
plt.imshow(LBT.wfs.cam.frame)
plt.colorbar()
plt.title('WFS Camera')
plt.show()

# Access components of the reference misregistration
print('Reference Shifts: X = ' + str(LBT.m_ref.shiftX) + ' m, Y = ' + str(LBT.m_ref.shiftY) + ' m')

#%% Initialize SPRINT
# This is required before you can run SPRINT, it determines the sensitivity matrices you'll use
# LBT analyser is hardcoded to just run SPRINT on shiftX, shiftY, and rotation 

LBT.init_SPRINT(mode=30,                      # 30 is the default mode but can be changed if desired
                recompute_sensitivity=True)   # if you've already generated sensitivity matrices for this configuration you may not need to do it again

#%% Define the raw LBT files you want to run SPRINT on, three files are needed to generate a demodulated signal with all the associated info

# Set the tracking number
trck = '20201001_075153'
# Set the folder where the files are located
loc = '../demodulated_slopes_orig/bin_'+str(LBT.binning)+'/'+str(trck)
# Save the full file path to each required file
slopes = loc + '/demodulated_slopes.fits'
phi = loc + '/phi.fits'
info = loc + '/data_info.fits'

#%% De-modulate the LBT data
LBT.get_on_sky_modulated_signal(slopes=slopes, phi=phi, info=info)

# This saves the 2D signal as LBT.slopes_2D, the 1D signal for SPRINT as LBT.on_sky_slopes, and the data info as LBT.data_info

#%% Run SPRINT on the demodulated signal and print the final estimates

# Plot the signal being sent to SPRINT in 2D
plt.figure()
plt.imshow(LBT.slopes_2D)
plt.title('Signal')
plt.show()

# Run SPRINT - SPRINT will run on whatever is saved as LBT.on_sky_slopes
LBT.run_SPRINT(n_iteration=6,         # set the number of iterations you want SPRINT to do, defaults to 3
               n_update_zero_point=0, # state the number of times you want to re-calculate your sensitivity matrices and update your zero-point. Default is 0
               precision=3,           # precision to round your estimates to, default is 3
               gain_estimation=1,     # gain to apply after one estimation
               dm_input=None)         # you can provide a dm_input, but this should always be None if you have param['isLBT'] = True


# Print the SPRINT estimates compared to the LBT flux estimates
print('Flux Method shift estimates = ' + str(round(LBT.data_info.sx,4)) + ',' + str(round(LBT.data_info.sy,4)))
print('SPRINT estimate (X shift, Y shift, rotation) = ' + str(LBT.misreg_est))
