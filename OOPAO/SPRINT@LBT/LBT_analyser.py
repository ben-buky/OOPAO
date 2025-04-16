# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 16:45:51 2025

@author: bbuky
"""

import matplotlib.pyplot as plt
import numpy as np
from OOPAO.tools.tools import emptyClass, read_fits
from OOPAO.tools.displayTools import displayMap
from OOPAO.MisRegistration import MisRegistration
from lbt_tools import get_int_mat_from_lbt
from build_LBT import build_LBT, ref_picker
from OOPAO.SPRINT import SPRINT

class LBT_analyser:
    
    def __init__(self,param,binning=1,misReg=None,psim=False,make_plots=True,n_modes=500,atm=True):
        
        """ LBT Analyser
        
        General class for all functions and variables required for analysis involving the LBT model, LBT data, and SPRINT.
        
        Initialization creates an LBT model using the build_LBT function, and an object containing the reference/zero-point mis-registration being used.
        
        """
        
        # store key variables within the class
        self.param = param
        self.binning = binning
        
        # set reference misreg for the system
        
        if misReg==None:
            self.m_ref = ref_picker(self.param,self.binning)
        else:
            self.m_ref = misReg
            
        self.tel, self.ngs, self.atm, self.dm_lbt, self.wfs, self.M2C_CL, self.calib_CL = build_LBT(self.param,self.binning,self.m_ref,psim,make_plots,n_modes,atm)
            
    def init_SPRINT(self,mode=30,n_mis_reg=3,recompute_sensitivity=True):
        
        # modal basis considered
        if np.isscalar(mode):
            index_modes = [mode]
        else:
            index_modes = mode
        basis =  emptyClass()
        basis.modes         = self.M2C_CL[:,index_modes]
        basis.extra         = 'LBT_KL_'+str(index_modes[0])              # EXTRA NAME TO DISTINGUISH DIFFERENT SENSITIVITY MATRICES, BE CAREFUL WITH THIS!     

        self.dm_lbt.coefs = basis.modes*1e-9
        self.tel.resetOPD()  
        self.ngs*self.tel*self.dm_lbt*self.wfs

        plt.figure()
        displayMap(self.tel.OPD)
        plt.title('KL Mode = ' + str(index_modes))
        plt.show()

        self.obj =  emptyClass()
        self.obj.ngs     = self.ngs
        self.obj.tel     = self.tel
        self.obj.atm     = self.atm
        self.obj.wfs     = self.wfs
        self.obj.dm      = self.dm_lbt
        self.obj.param   = self.param
            
        self.sprint = SPRINT(self.obj, basis, mis_registration_zero_point=self.m_ref, wfs_mis_registered=self.wfs, n_mis_reg=n_mis_reg, recompute_sensitivity=recompute_sensitivity)
        
    def get_on_sky_modulated_signal(self,slopes,phi,info):
        # raw slopes from LBT (scrambled)
        on_sky_slopes_raw = read_fits(slopes)
        phi_raw = read_fits(phi)

        # ordered slopes from LBT to comply with OOPAO logic
        on_sky_slopes_ordered = get_int_mat_from_lbt(self.param,self.wfs, IM_from_LBT = on_sky_slopes_raw)
        phi_ordered = get_int_mat_from_lbt(self.param,self.wfs, IM_from_LBT = phi_raw)
        
        # gather telemetry data from measurement
        data_info_raw = read_fits(info)

        # store the info into a Python class
        data_info = emptyClass()
        data_info.sx                = data_info_raw[0] # flux method estimate for X shift
        data_info.sy                = data_info_raw[1] # flux method estimate for Y shift
        data_info.re_rotator_angle  = data_info_raw[2]
        data_info.tel_rotator_angle = data_info_raw[3]
        data_info.mode_index        = int(data_info_raw[4])
        data_info.amplitude         = data_info_raw[5]
        data_info.frequency         = data_info_raw[6]
        data_info.n_iteration       = data_info_raw[7]
        
        reference_imat      = get_int_mat_from_lbt(self.param,self.wfs)
        reference_slope     = reference_imat[:,data_info.mode_index]

        # reconstruct the slopes from modulation parameter
        on_sky_slopes = (np.sin(phi_ordered))*on_sky_slopes_ordered/data_info.amplitude 
        n = len(on_sky_slopes)
        on_sky_slopes[:n//2]*=-1
        on_sky_slopes*= np.sign(reference_slope.T@on_sky_slopes/reference_slope.T@reference_slope)
        
        # create 2D slopes map using 2D map of valid pixels
        valid_pix_map = np.concatenate((self.wfs.userValidSignal,self.wfs.userValidSignal))
        valid_pix = np.reshape(valid_pix_map,valid_pix_map.shape[0]*valid_pix_map.shape[1]).astype('float64')
        # replace all the ones in the valid pixel array with the on_sky_slopes values
        valid_pix[valid_pix == 1] = on_sky_slopes
        self.slopes_2D = np.reshape(valid_pix,[valid_pix_map.shape[0],valid_pix_map.shape[1]])
        
        self.on_sky_slopes = on_sky_slopes
        self.data_info = data_info


    def run_SPRINT(self, n_iteration=3, n_update_zero_point=0, precision=3, gain_estimation=1, dm_input=None):
        
        self.sprint.estimate(self.obj, self.on_sky_slopes,n_iteration, n_update_zero_point, precision, gain_estimation, dm_input)
        
        # save the SPRINT estimates as a new variable
        self.misreg_est = self.sprint.mis_registration_buffer[-1]
        
    def init_SPRINT_ref(self,n_mis_reg=5,recompute_sensitivity=True):
        
        # modal basis considered
        index_modes = self.M2C_CL[:,:]
        basis =  emptyClass()
        basis.modes         = index_modes
        basis.extra         = 'LBT_KL_full'               # EXTRA NAME TO DISTINGUISH DIFFERENT SENSITIVITY MATRICES, BE CAREFUL WITH THIS!     

        #self.dm_lbt.coefs = basis.modes*1e-9
        #self.tel.resetOPD()  
        #self.tel*self.dm_lbt

        #plt.figure()
        #displayMap(self.tel.OPD)
        #plt.show()

        self.obj =  emptyClass()
        self.obj.ngs     = self.ngs
        self.obj.tel     = self.tel
        self.obj.atm     = self.atm
        self.obj.wfs     = self.wfs
        self.obj.dm      = self.dm_lbt
        self.obj.param   = self.param
            
        self.sprint = SPRINT(self.obj, basis, mis_registration_zero_point=self.m_ref, wfs_mis_registered=self.wfs, n_mis_reg=n_mis_reg, recompute_sensitivity=recompute_sensitivity)
        
    def run_ref(self,n_iteration=8, n_update_zero_point=0, precision=3, gain_estimation=1, dm_input=None):
        
        self.signal = self.calib_CL.D[:,:len(self.M2C_CL[1])]
        
        self.sprint.estimate(self.obj, self.signal, n_iteration, n_update_zero_point, precision, gain_estimation, dm_input)
        
        # save the SPRINT estimates
        self.ref_est = self.sprint.mis_registration_buffer[-1]
        