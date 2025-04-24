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

# #%
# zonal_im = LBT.calib_CL.D@np.linalg.pinv(LBT.M2C_CL)

# plt.close('all')
# ind = [0,1,2,3]
# ind =np.arange(0,100,4)

# zer = np.zeros(LBT.dm_lbt.nValidAct)
# zer[ind] = 1e-9

# display_wfs_signals(wfs=LBT.wfs, signals = np.sum(zonal_im[:,ind],axis=1))
# LBT.dm_lbt.coefs = zer
# LBT.tel*LBT.dm_lbt*LBT.wfs

# display_wfs_signals(wfs=LBT.wfs, signals = LBT.wfs.signal)

#%%
plt.close('all')
ind = np.arange(9)

LBT.dm_lbt.coefs = LBT.M2C_CL[:,ind]*1e-9
LBT.ngs*LBT.tel*LBT.dm_lbt*LBT.wfs


displayMap(LBT.tel.OPD)

a = display_wfs_signals(LBT.wfs, signals = LBT.calib_CL.D[:,ind],norma = True,returnOutput=True)
b = display_wfs_signals(LBT.wfs, signals = LBT.wfs.signal,norma=True,returnOutput=True) # calib_0 comes from computing meta sensitivity matrices

a[np.isinf(a)] =0
b[np.isinf(b)] =0

from OOPAO.tools.displayTools import interactive_show
plt.close('all')
interactive_show(a,b)

#%%
# pupil = LBT.tel.pupil.copy()
# #%%

# plt.close('all')

# ind = [5,6,7,8]
# # ind = [30]
# # ind = np.arange(49)
# ref_wfs=[]
# for i in range(4):
#     ref_wfs.append(display_wfs_signals(LBT.wfs, signals = LBT.calib_CL.D[:,ind[i]],norma = True, returnOutput=True))
# plt.close('all')
# for i in range(4):
#     LBT.dm_lbt.coefs = LBT.M2C_CL[:,ind[i]]*1e-9
#     LBT.ngs*LBT.tel*LBT.dm_lbt
#     ref_OPD_0 = LBT.tel.OPD.copy()
#     plt.figure()
#     #normal 
#     LBT.tel.pupil = pupil.copy()
#     ref_OPD = ref_OPD_0.copy()
#     plt.subplot(3,5,5)
#     plt.imshow(ref_wfs[i])
#     plt.subplot(3,5,15)
#     plt.imshow(ref_wfs[i])
#     plt.subplot(3,5,10)
#     plt.imshow(ref_wfs[i])

#     for i_r in range(4):
#         if i_r>0:
#             ref_OPD = np.rot90(ref_OPD) 
#             LBT.tel.pupil = np.rot90(LBT.tel.pupil)
#         LBT.tel.OPD = ref_OPD           
#         LBT.tel*LBT.wfs
#         signal_normal = LBT.wfs.signal.copy()
#         plt.subplot(3,5,i_r+1)
#         plt.imshow(LBT.wfs.signal_2D)
    
#     #flip 
#     ref_OPD = np.flip(ref_OPD_0.copy())
#     LBT.tel.pupil = np.flip(pupil.copy())

#     for i_r in range(4):
#         if i_r>0:
#             ref_OPD = np.rot90(ref_OPD)            
#             LBT.tel.pupil = np.rot90(LBT.tel.pupil)

#         LBT.tel.OPD = ref_OPD
#         LBT.tel*LBT.wfs
#         signal_normal = LBT.wfs.signal.copy()
#         plt.subplot(3,5,5+i_r+1)
#         plt.imshow(LBT.wfs.signal_2D)
        
#     #fliplr
#     ref_OPD = np.fliplr(ref_OPD_0.copy())
#     LBT.tel.pupil= np.fliplr(pupil.copy())

#     for i_r in range(4):
#         if i_r>0:
#             ref_OPD = np.rot90(ref_OPD)    
#             LBT.tel.pupil = np.rot90(LBT.tel.pupil)

#         LBT.tel.OPD = ref_OPD
#         LBT.tel*LBT.wfs
#         signal_normal = LBT.wfs.signal.copy()
#         plt.subplot(3,5,10+i_r+1)
#         plt.imshow(LBT.wfs.signal_2D)

        

# #     #normal 
# #     # LBT.tel*LBT.wfs
# #     # signal_normal = wfs.signal.copy()

    
    
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

#%%
from OOPAO.tools.tools import emptyClass
from OOPAO.SPRINT import SPRINT


# modal basis considered
index_modes = [30]
basis =  emptyClass()
basis.modes         = LBT.M2C_CL[:,index_modes]
basis.extra         = 'LBT_KL_'+str(index_modes[0])              # EXTRA NAME TO DISTINGUISH DIFFERENT SENSITIVITY MATRICES, BE CAREFUL WITH THIS!     
obj =  emptyClass()
obj.ngs     = LBT.ngs
obj.tel     = LBT.tel
obj.atm     = LBT.atm
obj.wfs     = LBT.wfs
obj.dm      = LBT.dm_lbt
obj.param   = LBT.param


sprint_dict = dict()

input_rot = [-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5]
input_rot = [-1.5,-0.5,0, 0.5,1.5]

for i in range(len(input_rot)):
    m = MisRegistration()
    m.rotationAngle = 299.515-0.5
    # m.rotationAngle = 299.515-0.5-5

    m.shiftX = 0.135
    m.shiftY = 0.005
    m.radialScaling = 0.021
    m.tangentialScaling = 0.026
    
    m.rotationAngle += input_rot[i]
    sprint_dict['sprint_'+str(i)] = SPRINT(obj,
                    basis,
                    mis_registration_zero_point=m,
                    wfs_mis_registered=LBT.wfs,
                    n_mis_reg=1,
                    recompute_sensitivity=True,
                    fast_algorithm = False,
                    ind_mis_reg=[2,0,1,3,4])



#%%

est = []
for i in range(len(input_rot)):
    plt.close('all')
    tmp_sprint = sprint_dict['sprint_'+str(i)]
    tmp_sprint.fast_algorithm = True
    tmp_sprint.estimate(obj,
                        on_sky_slopes = LBT.calib_CL.D[:,index_modes],
                        n_iteration=4,
                        tolerance=10)
    
    est.append(tmp_sprint.mis_registration_zero_point.rotationAngle -tmp_sprint.mis_registration_out.rotationAngle)
#%%
plt.figure(),
plt.plot(input_rot,input_rot,'-')
# plt.plot(input_rot,np.asarray(est).T,'-o'),
plt.plot(input_rot,np.asarray(est)-11.25,'-o'),

plt.xlabel('Input Rotation'),plt.ylabel('Output Rotation')
#%%
plt.figure(),
# plt.plot(input_rot,input_rot*0,'-')
plt.plot(input_rot,np.asarray(est) - np.asarray(input_rot),'-o'),
plt.xlabel('Input Rotation'),plt.ylabel('Output Rotation')
#%%
from OOPAO.tools.interpolateGeometricalTransformation import (anamorphosisImageMatrix, rotateImageMatrix, translationImageMatrix)
import skimage.transform as sk

def apply_mis_reg(tel,map_2d, misReg):
    pixelsize = tel.D/tel.resolution
    tel.resetOPD()
    # 2) transformations for the mis-registration
    anamMatrix              = anamorphosisImageMatrix(tel.OPD,misReg.anamorphosisAngle,[1+misReg.radialScaling,1+misReg.tangentialScaling])
    rotMatrix               = rotateImageMatrix(tel.OPD,misReg.rotationAngle)
    shiftMatrix             = translationImageMatrix(tel.OPD,[misReg.shiftY/pixelsize,misReg.shiftX/pixelsize]) #units are in m
          
    # 3) Global transformation matrix
    transformationMatrix    =  anamMatrix + rotMatrix + shiftMatrix 
    
    def globalTransformation(image):
            output  = sk.warp(image,(transformationMatrix).inverse,order=3)
            return output
    out = globalTransformation(map_2d)
    return out



LBT.dm_lbt.coefs = LBT.M2C_CL[:,30]
# LBT_.dm_lbt.coefs = LBT.M2C_CL[:,0]

OPD_in = LBT.tel.pupil.astype(float) #LBT.dm_lbt.OPD
# OPD_exp = LBT_.dm_lbt.OPD

misReg = MisRegistration()
misReg.rotationAngle = 90
for i in range(1):
    OPD_in = apply_mis_reg(LBT.tel, OPD_in, misReg)


interactive_show(LBT.tel.pupil.astype(float) ,OPD_in*LBT.tel.pupil)