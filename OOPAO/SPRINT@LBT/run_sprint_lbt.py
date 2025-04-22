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
from OOPAO.tools.displayTools import display_wfs_signals, displayMap
from lbt_tools import compare_wfs_signals
#%% read parameter file
from parameterFile_SOUL_I_Band_final import initializeParameterFile
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
# param['filename_m2c']       = loc + side+'/'+KL+'/m2c.fits'    

param['filename_coord']     = 'C:/Users/cheritier/OOPAO/user/SOUL/lbt_data/DATAfromLBT/left/KL_v20/act_coordinates.fits'

#%
from OOPAO.MisRegistration import MisRegistration


m = MisRegistration()
m.rotationAngle = 299.515 -0.5
m.shiftX = 0.135
m.shiftY = 0.005
m.radialScaling = 0.021
m.tangentialScaling = 0.026

# This initializes the class and creates the desired LBT model
LBT = LBT_analyser(param=param,     # define where to find the model parameters
                   binning=binning, # set your binning, must match what you use for BB_file_picker
                   misReg=m,     # you can provide a unique starting misregistration, None will use the standard zero-point reference for your binning
                   psim=False,      # you have the option to generate a synthetic IM rather than using the real one, the default is False
                   make_plots=True, # you can turn off the plots produced in the setup
                   n_modes=500,     # specify your number of modes, this will default to the maximum number allowed for your binning
                   atm=False)        # you can prevent an atmosphere being generated if you're only running SPRINT and want to speed up the code

# The above initialization stores your parameters as LBT.param and binning as LBT.binning. 
# You can also now access objects of the model and the reference mis-registration through your LBT, some examples are below:

# Display the telescope OPD 
# LBT.tel + LBT.atm # will only work if you've generated an atmosphere
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

#%%
plt.close('all')
ind = np.arange(9)

LBT.dm_lbt.coefs = LBT.M2C_CL[:,ind]*1e-9
LBT.ngs*LBT.tel*LBT.dm_lbt*LBT.wfs
compare_wfs_signals(LBT.wfs,LBT.calib_CL.D[:,ind],LBT.wfs.signal)
displayMap(LBT.tel.OPD)

#%% Initialize SPRINT
# This is required before you can run SPRINT, it determines the sensitivity matrices you'll use
# LBT analyser is hardcoded to just run SPRINT on shiftX, shiftY, and rotation 
# ind = [10,12,14,16,20,30,60,100,120]
ind = 30

LBT.init_SPRINT(mode=ind,                      # 30 is the default mode but can be changed if desired or be multiple modes
                n_mis_reg=3,                  # the number of mis-registration variables. The default is 3, which are shiftX, shiftY, and rotation
                recompute_sensitivity=True)   # if you've already generated sensitivity matrices for this configuration you may not need to do it again

# Example of using multiple modes for SPRINT:
#modes = [30,40,60,80]
#LBT.init_SPRINT(mode=modes,recompute_sensitivity=True)

#%%
LBT.on_sky_slopes = LBT.calib_CL.D[:,ind]
# Run SPRINT - SPRINT will run on whatever is saved as LBT.on_sky_slopes
LBT.run_SPRINT(n_iteration=5,         # set the number of iterations you want SPRINT to do, defaults to 3
               n_update_zero_point=0, # state the number of times you want to re-calculate your sensitivity matrices and update your zero-point. Default is 0
               precision=3,           # precision to round your estimates to, default is 3
               gain_estimation=1,     # gain to apply after one estimation
               dm_input=None,
               tolerance=10)         # you can provide a dm_input, but this should always be None if you have param['isLBT'] = True


#%%
compare_wfs_signals(LBT.wfs,LBT.calib_CL.D[:,ind],LBT.sprint.calib_last.D)

# compare_wfs_signals(LBT.wfs,LBT.calib_CL.D[:,ind],LBT.sprint.calib_0.D)
# display_wfs_signals(LBT.wfs, signals = LBT.calib_CL.D[:,ind],norma = True)
# display_wfs_signals(LBT.wfs, signals = LBT.sprint.calib_last.D,norma=True) # calib_0 comes from computing meta sensitivity matrices

#%% Initialize SPRINT
# This is required before you can run SPRINT, it determines the sensitivity matrices you'll use
# LBT analyser is hardcoded to just run SPRINT on shiftX, shiftY, and rotation 
# ind = [10,12,14,16,20,30,60,100,120]
ind = 30

LBT.init_SPRINT(mode=ind,                      # 30 is the default mode but can be changed if desired or be multiple modes
                n_mis_reg=3,                  # the number of mis-registration variables. The default is 3, which are shiftX, shiftY, and rotation
                recompute_sensitivity=True)   # if you've already generated sensitivity matrices for this configuration you may not need to do it again


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




# #%%
# from OOPAO.tools.tools import emptyClass
# from OOPAO.SPRINT import SPRINT


# # modal basis considered
# index_modes = [30]
# basis =  emptyClass()
# basis.modes         = LBT.M2C_CL[:,index_modes]
# basis.extra         = 'LBT_KL_'+str(index_modes[0])              # EXTRA NAME TO DISTINGUISH DIFFERENT SENSITIVITY MATRICES, BE CAREFUL WITH THIS!     
# obj =  emptyClass()
# obj.ngs     = LBT.ngs
# obj.tel     = LBT.tel
# obj.atm     = LBT.atm
# obj.wfs     = LBT.wfs
# obj.dm      = LBT.dm_lbt
# obj.param   = LBT.param


# sprint_dict = dict()

# input_rot = [-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5]
# input_rot = [-1.5,-0.5,0, 0.5,1.5]

# for i in range(len(input_rot)):
#     m = MisRegistration()
#     m.rotationAngle = 299.515-0.5
#     # m.rotationAngle = 299.515-0.5-5

#     m.shiftX = 0.135
#     m.shiftY = 0.005
#     m.radialScaling = 0.021
#     m.tangentialScaling = 0.026
    
#     m.rotationAngle += input_rot[i]
#     sprint_dict['sprint_'+str(i)] = SPRINT(obj,
#                     basis,
#                     mis_registration_zero_point=m,
#                     wfs_mis_registered=LBT.wfs,
#                     n_mis_reg=1,
#                     recompute_sensitivity=True,
#                     fast_algorithm = False,
#                     ind_mis_reg=[2,0,1,3,4])


# est = []
# for i in range(len(input_rot)):
#     plt.close('all')
#     tmp_sprint = sprint_dict['sprint_'+str(i)]
#     tmp_sprint.fast_algorithm = False
#     tmp_sprint.estimate(obj,
#                         on_sky_slopes = LBT.calib_CL.D[:,index_modes],
#                         n_iteration=4,
#                         tolerance=10)
    
#     est.append(tmp_sprint.mis_registration_zero_point.rotationAngle -tmp_sprint.mis_registration_out.rotationAngle)
    
# plt.figure(),
# plt.plot(input_rot,input_rot,'-')
# # plt.plot(input_rot,np.asarray(est).T,'-o'),
# plt.plot(input_rot,np.asarray(est)-11.25,'-o'),

# plt.xlabel('Input Rotation'),plt.ylabel('Output Rotation')