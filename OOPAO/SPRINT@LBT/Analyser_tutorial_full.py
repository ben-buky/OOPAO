# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 16:54:09 2025

@author: tsl29789

This tutorial aims to explain in as much detail as possible how to use LBT_analyser to run SPRINT on some LBT data.
"""
import matplotlib.pyplot as plt
import numpy as np
from LBT_analyser import LBT_analyser
from build_LBT import BB_file_picker
#%% read parameter file
from parameterFile_SOUL_I_Band import initializeParameterFile
param = initializeParameterFile()
#%% Create your LBT

# Set your binning, whether you're using new or old data, and which telescope you want 
# Then use BB_file_picker to generate the correct file paths parameters
param['new_IF'] = True
param['which_tel'] = 'Left'
binning=1
BB_file_picker(param, binning)

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

#%% Taissirs way to set the right parameters

param['new_IF'] = False 
side = 'left'
param['new_IF'] = True 
side = 'SX'

if param['new_IF']:
    loc = 'C:/Users/cheritier/OOPAO/user/SOUL/lbt_data/new_data_from_LBT/'
else:
    loc = 'C:/Users/cheritier/OOPAO/user/SOUL/lbt_data/old_data_from_LBT/'

    
if side == 'DX':
    trck_pup = '20200131_180431'
    trck_int_mat = '20221205_210843'
    KL = 'KL_v29'
    
if side == 'SX':
    trck_pup = '20250219_221329'
    trck_int_mat = '20181215_201757'
    KL = 'KL_v20'

if side == 'left':
    trck_pup = 'mode2/20190909_203854/'
    trck_int_mat = '20181215_201757'
    KL = 'KL_v20'    
    side = 'SX'

if side == 'right':
    trck_pup = '20250219_221329'
    trck_int_mat = '20181215_201757'
    KL = 'KL_v20'


param['filename_pup']       = loc + side+'/pupils/'+trck_pup+'/pup1.fits'
param['slopex']             = loc + side+'/pupils/'+trck_pup+'/slopex'
param['int_mat']            = loc + side+'/'+KL+'/RECs/Intmat_'+trck_int_mat+'.fits'
if param['new_IF']:
    param['filename_if']        = loc + side+'/'+KL+'/phase_matrix.sav'
    param['filename_mir_modes'] = loc + side+'/'+KL+'/phase_matrix.sav'
    param['filename_m2c']       = loc + side+'/'+KL+'/phase_matrix.sav'
else:
    param['filename_if']        = loc + side+'/'+KL+'/LBT672bIF.fits'
    param['filename_mir_modes'] = loc + side+'/'+KL+'/mirmodes.fits'
    param['filename_m2c']       = loc + side+'/'+KL+'/m2c.fits'    

param['filename_coord']     = 'C:/Users/cheritier/OOPAO/user/SOUL/lbt_data/DATAfromLBT/left/KL_v20/act_coordinates.fits'

#%%

# This initializes the class and creates the desired LBT model
LBT = LBT_analyser(param=param,     # define where to find the model parameters
                   binning=binning, # set your binning, must match what you use for BB_file_picker
                   misReg=None,     # you can provide a unique starting misregistration, None will use the standard zero-point reference for your binning
                   psim=False,      # you have the option to generate a synthetic IM rather than using the real one, the default is False
                   make_plots=True, # you can turn off the plots produced in the setup
                   n_modes=500,     # specify your number of modes, this will default to the maximum number allowed for your binning
                   atm=False)        # you can prevent an atmosphere being generated if you're only running SPRINT and want to speed up the code

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

LBT.init_SPRINT(mode=30,                      # 30 is the default mode but can be changed if desired or be multiple modes
                n_mis_reg=3,                  # the number of mis-registration variables. The default is 3, which are shiftX, shiftY, and rotation
                recompute_sensitivity=True)   # if you've already generated sensitivity matrices for this configuration you may not need to do it again

# Example of using multiple modes for SPRINT:
#modes = [30,40,60,80]
#LBT.init_SPRINT(mode=modes,recompute_sensitivity=True)

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

# Compare the true IM for mode 30 to what SPRINT has
from OOPAO.tools.displayTools import display_wfs_signals

display_wfs_signals(LBT.wfs, signals = LBT.calib_CL.D[:,30])
display_wfs_signals(LBT.wfs, signals = LBT.sprint.calib_0.D) # calib_0 comes from computing meta sensitivity matrices

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
    
    # Save the estimates
    Flux_ests[i] = [LBT.data_info.sx,LBT.data_info.sy]
    SPRINT_ests[i] = LBT.misreg_est

# Plot the results for X and Y shift
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
    
