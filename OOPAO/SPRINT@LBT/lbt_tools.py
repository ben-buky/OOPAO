# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:56:44 2021

@author: cheritie
"""
from astropy.io import fits as pfits
from joblib import Parallel, delayed
import scipy
import time
import numpy as np
import skimage.transform as sk
import matplotlib.pyplot as plt

from OOPAO.tools.interpolateGeometricalTransformation import rotateImageMatrix,rotation,translationImageMatrix,translation,anamorphosis,anamorphosisImageMatrix
from OOPAO.MisRegistration  import MisRegistration
from OOPAO.tools.tools import read_fits, emptyClass
from OOPAO.tools.displayTools import interactive_show,display_wfs_signals

def get_influence_functions_new(telescope, misReg, filename_IF,filename_mir_modes,filename_coordinates, filename_M2C):
    # misReg.rotationAngle += 0
    #% Load eigen modes of the mirror
    tmp = scipy.io.readsav(filename_IF)
    influenceFunctions_ASM = np.zeros([tmp['dpix']*tmp['dpix'],tmp['klmatrix'].shape[1]])
    influenceFunctions_ASM[tmp['idx_mask'],:] = tmp['klmatrix']
    influenceFunctions_ASM = influenceFunctions_ASM.reshape(tmp['dpix'],tmp['dpix'],influenceFunctions_ASM.shape[1]).T
    # influenceFunctions_ASM = np.moveaxis(influenceFunctions_ASM,1,2)

    # hdu = pfits.open(filename_IF)
    # F = hdu[0].data.byteswap().newbyteorder()
    # influenceFunctions_ASM = F[:,30:470,19:459]
    a= time.time()
    nAct,nx, ny = influenceFunctions_ASM.shape
    pixelSize_ASM_original = 8.25/nx
    # size of the influence functions maps
    resolution_ASM_original       = int(nx)   
    # resolution of the M1 pupil
    resolution_M1                = telescope.pupil.shape[1]
    # compute the pixel scale of the M1 pupil
    pixelSize_M1                 = 8.25/resolution_M1
    # compute the ratio_ASM_M1 between both pixel scale.
    ratio_ASM_M1                  = pixelSize_ASM_original/pixelSize_M1
    # after the interpolation the image will be shifted of a fraction of pixel extra if ratio_ASM_M1 is not an integer
    extra = (ratio_ASM_M1)%1 

    # difference in pixels between both resolutions    
    nPix = resolution_ASM_original-resolution_M1
    
    if nPix%2==0:
        # case nPix is even
        # alignement of the array with respect to the interpolation 
        # (ratio_ASM_M1 is not always an integer of pixel)
        extra_x = extra/2 -0.5
        extra_y = extra/2 -0.5
        
        # crop one extra pixel on one side
        nCrop_x = nPix//2
        nCrop_y = nPix//2
    else:
        # case nPix is uneven
        # alignement of the array with respect to the interpolation 
        # (ratio_ASM_M1 is not always an integer of pixel)
        extra_x = extra/2 -0.5 -0.5
        extra_y = extra/2 -0.5 -0.5
        # crop one extra pixel on one side
        nCrop_x = nPix//2
        nCrop_y = (nPix//2)+1
           
    # allocate memory to store the influence functions
    influMap = np.zeros([resolution_ASM_original,resolution_ASM_original])  
    
    #-------------------- The Following Transformations are applied in the following order -----------------------------------
       
    # 1) Down scaling to get the right pixel size according to the resolution of M1
    downScaling     = anamorphosisImageMatrix(influMap,0,[ratio_ASM_M1,ratio_ASM_M1])
    
    # 2) transformations for the mis-registration
    anamMatrix              = anamorphosisImageMatrix(influMap,misReg.anamorphosisAngle,[1+misReg.radialScaling,1+misReg.tangentialScaling])
    rotMatrix               = rotateImageMatrix(influMap,misReg.rotationAngle)
    shiftMatrix             = translationImageMatrix(influMap,[-misReg.shiftX/pixelSize_M1,-misReg.shiftY/pixelSize_M1]) #units are in m
    
    # Shift of half a pixel to center the images on an even number of pixels
    alignmentMatrix         = translationImageMatrix(influMap,[extra_x,extra_y])
        
    # 3) Global transformation matrix
    transformationMatrix    = downScaling + anamMatrix + rotMatrix + shiftMatrix + alignmentMatrix
    
    def globalTransformation(image):
            output  = sk.warp(image,(transformationMatrix).inverse,order=0)
            return output
    # definition of the function that is run in parallel for each 
    def reconstruction_IF(influMap):
        output = globalTransformation(np.fliplr(influMap))  
        # output = globalTransformation(influMap)

        return output    
    def joblib_reconstruction():
        Q=Parallel(n_jobs=4,prefer='threads')(delayed(reconstruction_IF)(i) for i in influenceFunctions_ASM)
        return Q 

    influenceFunctions_tmp =  np.moveaxis(np.asarray(joblib_reconstruction()),0,-1)
    influenceFunctions = influenceFunctions_tmp  [nCrop_x:-nCrop_y,nCrop_x:-nCrop_y,:]
    influenceFunctions = -influenceFunctions.reshape(influenceFunctions.shape[0]*influenceFunctions.shape[1],nAct)
    
    # b= time.time()
    # hdu = pfits.open(filename_mir_modes)
    # U = hdu[0].data
    # mod2zon = np.reshape(U[np.where(np.abs(U)!=0)],[663,663])
    # influenceFunctions = influenceFunctions_tmp@ mod2zon.T
    # if filename_IF == filename_M2C:
    #     M2C = tmp['klm2c']
    #     validAct = np.where(M2C[:,2]!=0)
    #     validAct = validAct[0].astype(int)
    #     M2C = M2C[validAct,:influenceFunctions.shape[2]]
    #     M2C = np.eye(nAct)
    # else:
    M2C = read_fits(filename_M2C)
    validAct = np.where(M2C[:,2]!=0)
    validAct = validAct[0].astype(int)
    M2C = M2C[validAct,:]
    # M2C *=-1

    # print(M2C.shape)
    # print(influenceFunctions.shape)
    influenceFunctions = influenceFunctions@np.linalg.pinv(M2C[:,:influenceFunctions.shape[1]])
    
    influenceFunctions = influenceFunctions.reshape(resolution_M1,resolution_M1,M2C.shape[0])

    hdu = pfits.open(filename_coordinates)
    coordinates_ASM_original =   hdu[0].data[validAct,:]/100
    # recenter the initial coordinates_M4_originalinates around 0
    coordinates_ASM = (coordinates_ASM_original)*ratio_ASM_M1  
    return influenceFunctions, coordinates_ASM, M2C, validAct


def get_influence_functions(telescope, misReg, filename_IF,filename_mir_modes,filename_coordinates, filename_M2C):
    #% Load eigen modes of the mirror
    hdu = pfits.open(filename_IF)
    F = hdu[0].data.byteswap().newbyteorder()
    influenceFunctions_ASM = F[:,30:470,19:459]
    # influenceFunctions_ASM = np.moveaxis(influenceFunctions_ASM, 1,2)
    a= time.time()
    nAct,nx, ny = influenceFunctions_ASM.shape
    pixelSize_ASM_original = 8.25/nx
    # size of the influence functions maps
    resolution_ASM_original       = int(nx)   
    # resolution of the M1 pupil
    resolution_M1                = telescope.pupil.shape[1]
    # compute the pixel scale of the M1 pupil
    pixelSize_M1                 = 8.25/resolution_M1
    # compute the ratio_ASM_M1 between both pixel scale.
    ratio_ASM_M1                  = pixelSize_ASM_original/pixelSize_M1
    # after the interpolation the image will be shifted of a fraction of pixel extra if ratio_ASM_M1 is not an integer
    extra = (ratio_ASM_M1)%1 

    # difference in pixels between both resolutions    
    nPix = resolution_ASM_original-resolution_M1
    
    if nPix%2==0:
        # case nPix is even
        # alignement of the array with respect to the interpolation 
        # (ratio_ASM_M1 is not always an integer of pixel)
        extra_x = extra/2 -0.5
        extra_y = extra/2 -0.5
        
        # crop one extra pixel on one side
        nCrop_x = nPix//2
        nCrop_y = nPix//2
    else:
        # case nPix is uneven
        # alignement of the array with respect to the interpolation 
        # (ratio_ASM_M1 is not always an integer of pixel)
        extra_x = extra/2 -0.5 -0.5
        extra_y = extra/2 -0.5 -0.5
        # crop one extra pixel on one side
        nCrop_x = nPix//2
        nCrop_y = (nPix//2)+1
           
    # allocate memory to store the influence functions
    influMap = np.zeros([resolution_ASM_original,resolution_ASM_original])  
    
    #-------------------- The Following Transformations are applied in the following order -----------------------------------
       
    # 1) Down scaling to get the right pixel size according to the resolution of M1
    downScaling     = anamorphosisImageMatrix(influMap,0,[ratio_ASM_M1,ratio_ASM_M1])
    
    # 2) transformations for the mis-registration
    anamMatrix              = anamorphosisImageMatrix(influMap,misReg.anamorphosisAngle,[1+misReg.radialScaling,1+misReg.tangentialScaling])
    rotMatrix               = rotateImageMatrix(influMap,misReg.rotationAngle)
    shiftMatrix             = translationImageMatrix(influMap,[-misReg.shiftX/pixelSize_M1,-misReg.shiftY/pixelSize_M1]) #units are in m
    
    # Shift of half a pixel to center the images on an even number of pixels
    alignmentMatrix         = translationImageMatrix(influMap,[extra_x,extra_y])
        
    # 3) Global transformation matrix
    transformationMatrix    = downScaling + anamMatrix + rotMatrix + shiftMatrix + alignmentMatrix
    
    def globalTransformation(image):
            output  = sk.warp(image,(transformationMatrix).inverse,order=0)
            return output
    # definition of the function that is run in parallel for each 
    def reconstruction_IF(influMap):
        output = globalTransformation(influMap)  
        return output    
    def joblib_reconstruction():
        Q=Parallel(n_jobs=4,prefer='threads')(delayed(reconstruction_IF)(i) for i in influenceFunctions_ASM)
        return Q 
    influenceFunctions_tmp =  np.moveaxis(np.asarray(joblib_reconstruction()),0,-1)
    influenceFunctions_tmp = influenceFunctions_tmp  [nCrop_x:-nCrop_y,nCrop_x:-nCrop_y,:]
    b= time.time()
    M2C = read_fits(filename_M2C)
    validAct = np.where(M2C[:,2]!=0)
    validAct = validAct[0].astype(int)
    M2C = M2C[validAct,:]
    
    hdu = pfits.open(filename_mir_modes)
    U = hdu[0].data
    mod2zon = np.reshape(U[np.where(np.abs(U)!=0)],[influenceFunctions_tmp.shape[2],influenceFunctions_tmp.shape[2]])
    influenceFunctions = -influenceFunctions_tmp@ mod2zon.T

    # coordinates of ASM before the interpolation        
    hdu = pfits.open(filename_coordinates)
    coordinates_ASM_original =   hdu[0].data[validAct,:]/100
    # recenter the initial coordinates_M4_originalinates around 0
    coordinates_ASM = (coordinates_ASM_original)*ratio_ASM_M1    
    
    # apply the transformations and re-center them for the new resolution resolution_M1
    coordinates_M4 = translation(rotation(anamorphosis(coordinates_ASM,misReg.anamorphosisAngle*np.pi/180,misReg.radialScaling,misReg.tangentialScaling),misReg.rotationAngle*np.pi/180),[misReg.shiftX/pixelSize_M1,misReg.shiftY/pixelSize_M1])+resolution_M1/2
    
    
    return influenceFunctions, coordinates_ASM_original, M2C, validAct


def get_wfs_pupil(filename,N0 = 120,binning =1):
    N1 = int(N0)
    N2 = int(N1*2/binning) 
    A= np.zeros(N2*N1)
    hdu = pfits.open(filename)
    index_pup = hdu[0].data
    index_pup=index_pup.astype(int)
    A[index_pup]=1
    wfsPup_ = np.reshape(A,[N1,N2])
    lin,col = np.where(wfsPup_>0)
    wfsPup=wfsPup_[min(lin):max(lin)+1,min(col):max(col)+1]
    wfsPup = np.flip(np.fliplr(np.rot90(((wfsPup)))))
    
    return wfsPup


def get_int_mat_from_lbt(param, wfs, IM_from_LBT = None):
    wfs_pup = wfs.userValidSignal.astype(float)
    binning = wfs.binning
    
    hdu = pfits.open(param['int_mat'])
    filename = param['slopex']

    with open(filename) as f:
        mylist = (f.read().splitlines())
    slopes_frame = (np.asarray(mylist).astype(int))
    #test1 = np.where(slopes_frame!=-1)
    #print(test1[0].shape)
    
    #!!! ordering of the slopes is [sx,sy,sx,sy,...]
    pix2slopes = [np.where(slopes_frame!=-1)]
    if IM_from_LBT is None:
        IM_from_LBT = hdu[0].data
        print(IM_from_LBT.shape)
        
        
    if np.ndim(IM_from_LBT)==1: 
        valid_pix       = (np.where(abs(IM_from_LBT)>0))            
        IM_for_model_0  = np.squeeze(IM_from_LBT[valid_pix])
        N1 = 120
        N2 = int(N1*2/binning)
        
        A= np.zeros([N2*N1])
        A[pix2slopes]= IM_for_model_0[slopes_frame[pix2slopes]]
    
        valid_signals_T = np.concatenate((wfs_pup,wfs_pup)).T    
        valid_signals = np.concatenate((wfs_pup,wfs_pup))
        IM_for_model = IM_for_model_0.copy()
            
        tmp = np.copy(valid_signals_T)
        tmp[np.where(valid_signals_T>0)] = A[pix2slopes]
        tmp_T = tmp.T
        IM_for_model[:] = tmp_T[np.where(valid_signals>0)]
        nSlopes = IM_for_model_0.shape[0]
        IM_for_model[1+nSlopes//2 : -1] = -IM_for_model[1+nSlopes//2 : -1] 

    else:
        valid_pix       = (np.where(abs(IM_from_LBT[:,10])>0))    
        IM_for_model_0  = np.squeeze(IM_from_LBT[valid_pix,:])
        N1 = 120
        N2 = int(N1*2/binning)
        
        A= np.zeros([N2*N1,IM_for_model_0.shape[1]])
    
        A[pix2slopes,:]= IM_for_model_0[slopes_frame[pix2slopes],:]
    
       
        valid_signals_T = np.concatenate((wfs_pup,wfs_pup)).T
        #print(valid_signals_T)
        
        valid_signals = np.concatenate((wfs_pup,wfs_pup))
    
        IM_for_model = IM_for_model_0.copy()
            
        for i in range(IM_for_model_0.shape[1]):
            tmp = np.copy(valid_signals_T)
           # test = tmp[np.where(valid_signals_T>0)]
            #print(test.shape)
            tmp[np.where(valid_signals_T>0)] = A[pix2slopes,i]
            tmp_T = tmp.T
            IM_for_model[:,i] = tmp_T[np.where(valid_signals>0)]
            
        nSlopes = IM_for_model_0.shape[0]
        
        IM_for_model[1+nSlopes//2 : -1,:] = -IM_for_model[1+nSlopes//2 : -1,:]

    return IM_for_model




def get_on_sky_modulated_signal(param,wfs,trck):
    # raw slopes from LBT (scrambled)
    on_sky_slopes_raw = read_fits('C:/Users/tsl29789/OneDrive - Science and Technology Facilities Council/Documents/AO/LBT/SPRINT@LBT/lbt_data/demodulated_slopes/bin_'+str(wfs.binning)+'/'+str(trck)+'/demodulated_slopes.fits')
    phi_raw = read_fits('C:/Users/tsl29789/OneDrive - Science and Technology Facilities Council/Documents/AO/LBT/SPRINT@LBT/lbt_data/demodulated_slopes/bin_'+str(wfs.binning)+'/'+str(trck)+'/phi.fits')

    # ordered slopes from LBT to comply with OOPAO logic
    on_sky_slopes_ordered = get_int_mat_from_lbt(param,wfs, IM_from_LBT = on_sky_slopes_raw)
    phi_ordered = get_int_mat_from_lbt(param,wfs, IM_from_LBT = phi_raw)
    
    # gather telemetry data from measurement
    data_info_raw = read_fits('C:/Users/tsl29789/OneDrive - Science and Technology Facilities Council/Documents/AO/LBT/SPRINT@LBT/lbt_data/demodulated_slopes/bin_'+str(wfs.binning)+'/'+str(trck)+'/data_info.fits')

    # store the info into a Python class
    data_info = emptyClass()
    data_info.sx                = data_info_raw[0]
    data_info.sy                = data_info_raw[1]
    data_info.re_rotator_angle  = data_info_raw[2]
    data_info.tel_rotator_angle = data_info_raw[3]
    data_info.mode_index        = int(data_info_raw[4])
    data_info.amplitude         = data_info_raw[5]
    data_info.frequency         = data_info_raw[6]
    data_info.n_iteration       = data_info_raw[7]
    
    reference_imat      = get_int_mat_from_lbt(param,wfs)
    reference_slope     = reference_imat[:,data_info.mode_index]

    # reconstruct the slopes from modulation parameter
    on_sky_slopes = (np.sin(phi_ordered))*on_sky_slopes_ordered/data_info.amplitude 
    n = len(on_sky_slopes)
    on_sky_slopes[:n//2]*=-1
    on_sky_slopes*= np.sign(reference_slope.T@on_sky_slopes/reference_slope.T@reference_slope)
    
    # create 2D slopes map using 2D map of valid pixels
    valid_pix_map = np.concatenate((wfs.userValidSignal,wfs.userValidSignal))
    valid_pix = np.reshape(valid_pix_map,valid_pix_map.shape[0]*valid_pix_map.shape[1]).astype('float64')
    # replace all the ones in the valid pixel array with the on_sky_slopes values
    valid_pix[valid_pix == 1] = on_sky_slopes
    slopes_2D = np.reshape(valid_pix,[valid_pix_map.shape[0],valid_pix_map.shape[1]])

    return on_sky_slopes, slopes_2D, data_info

def compare_wfs_signals(wfs,signal_1,signal_2,norma=True):
    a = display_wfs_signals(wfs, signals = signal_1, norma = norma,returnOutput=True)
    b = display_wfs_signals(wfs, signals = signal_2, norma = norma,returnOutput=True) # calib_0 comes from computing meta sensitivity matrices

    a[np.isinf(a)] =0
    b[np.isinf(b)] =0

    plt.close('all')
    interactive_show(a,b)
    
    return