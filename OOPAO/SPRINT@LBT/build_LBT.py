# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:06:11 2025

@author: bbuky
"""

import matplotlib.pyplot as plt
import numpy as np
#import time
#from astropy.io import fits as pfits
from OOPAO.MisRegistration import MisRegistration
from OOPAO.Atmosphere import Atmosphere
from OOPAO.Pyramid import Pyramid
#from OOPAO.DeformableMirror import DeformableMirror
from OOPAO.Telescope import Telescope
from OOPAO.Source import Source
#from OOPAO.calibration.ao_calibration import ao_calibration
from OOPAO.calibration.CalibrationVault import CalibrationVault
from OOPAO.calibration.InteractionMatrix import InteractionMatrix
from OOPAO.mis_registration_identification_algorithm.applyMisRegistration import applyMisRegistration
#from joblib import Parallel, delayed
#import skimage.transform as sk
#import ctypes
from OOPAO.tools.displayTools import displayMap#, makeSquareAxes
from lbt_tools import get_wfs_pupil, get_int_mat_from_lbt

def BB_file_picker(param,which_tel,binning):
    
    # set whether you're using new or old data
    if param['new_IF']:
        loc = '../new_data_from_lbt/'
        
        # determine which telescope you're using
        if which_tel == 'Left':
            
            side = 'SX'
            KL = 'KL_v20'
            # tracking numbers for the valid pupil mask for the PWFS
            trck_pup = ['bin1/20250219_221329','bin2/20250222_202110','bin3/20250222_212311','bin4/20230320_000004']
            # tracking numbers for the interaction matrices - binning=3 doesn't work and binning=4 doesn't exist hence 444444
            trck_int_mat = ['20181215_201757','20190401_211243','20200225_144315','bin4_old']
            
        if which_tel == 'Right':
            
            side = 'DX'
            KL = 'KL_v29'
            # tracking numbers for the valid pupil mask for the PWFS
            trck_pup = ['bin1/20200131_180431','bin2/20200205_211349','bin3/20230117_192316','bin4/20230221_093800']
            # tracking numbers for the interaction matrices - binning=3 doesn't work
            trck_int_mat = ['20230117_165253','20230120_192846','20230117_224144','20230118_200801']
         
        # set parameters
        param['filename_pup']  = loc + side+'/pupils/'+trck_pup[binning-1]+'/pup1.fits'
        param['slopex']        = loc + side+'/pupils/'+trck_pup[binning-1]+'/slopex'
        param['int_mat']       = loc + side+'/'+KL+'/RECs/Intmat_'+trck_int_mat[binning-1]+'.fits'
        
        # influence functions file
        param['filename_if']        = loc + side+'/'+KL+'/phase_matrix.sav'
        # ASM eigen modes
        param['filename_mir_modes'] = loc + side+'/'+KL+'/phase_matrix.sav'
        # Mode to Command matrix
        param['filename_m2c']       = loc + side+'/'+KL+'/phase_matrix.sav'
        
        # all configurations use the same ASM coordinates file but I've made two copies of it
        param['filename_coord'] = loc + 'act_coordinates.fits'
    
    else:
        # the old data only covers the left telescope
        loc = '../old_data_from_lbt'
        
        # tracking numbers
        trck_pup = ['mode2/20190909_203854','mode3/20190606_074345','mode4/20200210_000000','mode5/20200916_173857']
        trck_int_mat = ['20181215_201757','20190401_211243','20200225_144315','20200917_152928']
        
        # set parameters
        param['filename_pup']  = loc + '/SX/pupils/'+trck_pup[binning-1]+'/pup1.fits'
        param['slopex']        = loc + '/SX/pupils/'+trck_pup[binning-1]+'/slopex.dat'
        param['int_mat']       = loc + '/SX/KL_v20/RECs/Intmat_'+trck_int_mat[binning-1]+'.fits'
        
        param['filename_if']        = loc + '/SX/KL_v20/LBT672bIF.fits'
        param['filename_mir_modes'] = loc + '/SX/KL_v20/mirmodes.fits'
        param['filename_m2c']       = loc + '/SX/KL_v20/m2c.fits'
        # all configurations use the same ASM coordinates file but I've made two copies of it
        param['filename_coord'] = loc + '/act_coordinates.fits'


def build_LBT(param,binning=1,misReg=None,psim=False,make_plots=True,n_modes=500,atm=True):
    
    """
    This function creates the full LBT model from the real influence functions and a parameter file.

    Parameters
    ----------
    param : dict
        A dictionary which specifies all of the key parameters for the model.
    binning : float, optional
        Binning size for wavefront sensor.
        The default is 1. 
    misReg : float, optional
        The reference misregistration to be applied to the model as its zeropoint.
        The default is None; this will lead to the standard reference calculated by SPRINT to be used.
    psim : bool, optional
        Gives the user the option to calibrate the system and use a pseudo synthetic IM rather than the real LBT IM.
        The default is False, meaning the real IM will be used.
    make_plots : bool, optional
        Produces plots for each element of the system (eg. atmosphere, telescope etc.) as they are created.
        The default is True.
    n_modes : float, optional
        Defines the number of KL modes usd in the system.
        The default is 500; but will ultimately match the size of the real IM for all binnings unless a different number is specified.
    atm: bool, optional
        Determines whether an atmosphere should be generated or not.
        The default is True.

    Returns
    -------
    tel : telescope object
    ngs : natural guide star object 
    atm : atmosphere object
    dm : deformable mirror object 
    wfs : wavefront sensor object
    M2C_KL : M2C matrix with shape nAct vs nModes
    calib_modal : CalibrationVault object containing the interaction and command matrices for the system

    """
    
    # -----------------------    TELESCOPE   ----------------------------------
    tel = Telescope(resolution = param['resolution'],
                    diameter = param['diameter'],
                    samplingTime = param['samplingTime'],
                    centralObstruction = param['centralObstruction'])
    
    if make_plots is True:
        plt.figure()
        plt.imshow(tel.pupil)
        plt.show()
        
    # -----------------------     SOURCE   ----------------------------------
    ngs=Source(optBand   = param['opticalBand'],\
               magnitude = param['magnitude'])

    # combine the NGS to the telescope using '*' operator:
    ngs*tel
    tel.computePSF(zeroPaddingFactor = 8)
    if make_plots is True:
        plt.figure()
        plt.imshow(np.log(np.abs(tel.PSF)/np.max(np.abs(tel.PSF))),extent = [tel.xPSF_arcsec[0],tel.xPSF_arcsec[1],tel.xPSF_arcsec[0],tel.xPSF_arcsec[1]])
        plt.clim([-12,0])
        plt.show()
        
    # -----------------------     ATMOSPHERE   ----------------------------------
    if atm is True:
        atm=Atmosphere(telescope     = tel,\
                       r0            = param['r0'],\
                       L0            = param['L0'],\
                       windSpeed     = param['windSpeed'],\
                       fractionalR0  = param['fractionnalR0'],\
                       windDirection = param['windDirection'],\
                       altitude      = param['altitude'])
        # initialize atmosphere
        atm.initializeAtmosphere(tel)
    else:
        atm=None 
        
    # --------------------     REFERENCE MISREG   -------------------------------
    
    # NEED SEPARATE MISREG OBJECT JUST FOR DM SO WE HAVE NO SHIFTS WHEN CALCULATING M2C
    
    if misReg is None:
        # apply standard reference misreg from full SPRINT calculation
        m_ref_wfs = MisRegistration()     
        m_ref_wfs.shiftX            = 0.137   # in metres
        m_ref_wfs.shiftY            = -0.004  # in metres
        
        m_ref_dm = MisRegistration()
        m_ref_dm.rotationAngle     = 299.486 # in degrees
        m_ref_dm.radialScaling     = 0.025
        m_ref_dm.tangentialScaling = 0.019
        
    else:
        m_ref_wfs = MisRegistration()     
        m_ref_wfs.shiftX            = misReg.shiftX
        m_ref_wfs.shiftY            = misReg.shiftY
        
        m_ref_dm = MisRegistration()
        m_ref_dm.rotationAngle     = misReg.rotationAngle
        m_ref_dm.radialScaling     = misReg.radialScaling
        m_ref_dm.tangentialScaling = misReg.tangentialScaling
    
    # ----------------   BUILD INFLUENCE FUNCTIONS   ----------------------------
    
    if param['new_IF']:
        from lbt_tools import get_influence_functions_new as get_influence_functions
    else:
        from lbt_tools import get_influence_functions as get_influence_functions
    
    # compute the LBT ASM influence function applying the proper sampling and mis-registrations
    modes_lbt, coord_lbt, M2C, validAct =  get_influence_functions(telescope = tel,
                                                                   misReg = m_ref_dm,
                                                                   filename_IF = param['filename_if'],
                                                                   filename_mir_modes = param['filename_mir_modes'],
                                                                   filename_coordinates = param['filename_coord'],
                                                                   filename_M2C = param['filename_m2c'])
    
    if param['new_IF']:
        dm_sum = np.reshape(np.sum(modes_lbt**2, axis =2),[tel.resolution, tel.resolution])
        tel.pupil = (dm_sum!=0)
    
    # ---------------------     GET PUPIL MASK   --------------------------------
    
    filename_pup = param['filename_pup']

    mask = get_wfs_pupil(filename_pup, binning = binning)

    if make_plots is True:
        plt.figure()
        plt.imshow(mask)
        plt.show()

    # -------------------     WAVEFRONT SENSOR   --------------------------------
    # create the Pyramid Object
    wfs = Pyramid(nSubap                = param['nSubaperture'],\
                  telescope             = tel,\
                  modulation            = param['modulation'],\
                  lightRatio            = 0.1,\
                  n_pix_separation      = 8,\
                  n_pix_edge            = 16,\
                  calibModulation       = 0,\
                  psfCentering          = False,\
                  postProcessing        = param['postProcessing'],
                  userValidSignal       = mask.astype(int), 
                  binning               = binning) 
    
    if make_plots is True:
        plt.figure()
        plt.imshow(wfs.referenceSignal_2D)
        plt.colorbar()
        plt.title('Reference Signal before shifting')
        plt.show()
    
        plt.figure()
        plt.imshow(wfs.cam.frame)
        plt.colorbar()
        print(wfs.cam.frame.shape)
        plt.title('Camera Frame before shifting')
        plt.show()
    
        plt.figure()
        plt.imshow(wfs.grabQuadrant(1))
        plt.colorbar()
        plt.title('Quadrant 1 before shifting')
        plt.show()
        
    
    
    # -------------------     DEFORMABLE MIRROR   -------------------------------
    # WARNING: wfs is provided to the applyMisRegistration, the WFS pupils are shifted instead of the DM. 
    # this is important for the LBT model because the ASM is the pupil of the system. 
    # WFS must exist before misreg applied
    dm = applyMisRegistration(tel = tel,
                              misRegistration_tmp = m_ref_wfs+m_ref_dm,
                              param = param,
                              wfs = wfs)

    if make_plots is True:
        
        # Update WFS plots after shifting
        plt.figure()
        plt.imshow(wfs.cam.frame)
        plt.colorbar()
        print(wfs.cam.frame.shape)
        plt.title('Camera Frame after shifting')
        plt.show()
    
        plt.figure()
        plt.imshow(wfs.grabQuadrant(1))
        plt.colorbar()
        plt.title('Quadrant 1 after shifting')
        plt.show()
        
        # apply the KL Modes on the DM to check the quality of the model
        tel.resetOPD()
        dm.coefs = M2C[:, :10]
        tel*dm
        displayMap(tel.OPD, norma = True)
        plt.show()

        # display quadratic sum of IFs to show the actuators positions wrt the pupil
    
        dm_sum = np.reshape(np.sum(dm.modes**2, axis =1),[tel.resolution, tel.resolution])
        dm_sum /=dm_sum.max()
    
        plt.figure()
        plt.imshow(dm_sum + 0.5*tel.pupil, extent = [-tel.D/2, tel.D/2,-tel.D/2, tel.D/2])
        plt.show()
    
    # -------------------    Interaction Matrix   -------------------------------
    # case where we want synthetic IM and need to do calibration
    if psim is True:
        calib_modal =  InteractionMatrix(ngs            = ngs,
                                        atm            = atm,
                                        tel            = tel,
                                        dm             = dm,
                                        wfs            = wfs,
                                        M2C            = M2C[:,:n_modes], # M2C matrix used
                                        stroke         = 1e-9,    # Stroke for the push/pull in M2C units
                                        nMeasurements  = 6,         # Number of simultaneous measurements
                                        noise          = 'off',     # Disable wfs.cam noise
                                        display        = True,      # Display the time using tqdm
                                        single_pass    = True)      # Only push to compute the interaction matrix instead of push-pull
    
    # case where we're using the real LBT IM
    else:
        int_mat = get_int_mat_from_lbt(param,wfs)
        
        # if using binning not equal to 1, let the interaction matrix define the number of modes, unless the user has already defined a suitably small number of modes
        if binning != 1 and n_modes > len(int_mat[1]):
            n_modes = len(int_mat[1])
        
        calib_modal    = CalibrationVault(int_mat[:,:n_modes])
        
    # crop M2C matrix to correct size
    M2C_KL = M2C[:,:n_modes]

    return  tel, ngs, atm, dm, wfs, M2C_KL, calib_modal

# def ref_picker(param,)