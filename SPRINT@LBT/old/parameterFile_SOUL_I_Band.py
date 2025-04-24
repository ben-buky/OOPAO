# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:36:02 2020

@author: cheritie

THIS FILE WAS USED WITH THE REAL IM FOR SPRINT REFERENCE CALCULATIONS
"""


from OOPAO.tools.tools  import createFolder

def initializeParameterFile():
    # initialize the dictionaries
    param = dict()
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATMOSPHERE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['r0'                   ] = 0.15                                           # value of r0 in the visibile in [m]
    param['L0'                   ] = 30                                             # value of L0 in the visibile in [m]
    param['fractionnalR0'        ] = [1]                                    # Cn2 profile
    param['windSpeed'            ] = [10]                                      # wind speed of the different layers in [m.s-1]
    param['windDirection'        ] = [120]                                     # wind direction of the different layers in [degrees]
    param['altitude'             ] = [1000]                                   # altitude of the different layers in [m]
    
                              
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M1 PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['diameter'             ] = 8.25                                             # diameter in [m]
    param['nSubaperture'         ] = 40                                             # number of PWFS subaperture along the telescope diameter
    param['nPixelPerSubap'       ] = 4                                            # sampling of the PWFS subapertures
    param['resolution'           ] = param['nSubaperture']*param['nPixelPerSubap']  # resolution of the telescope driven by the PWFS
    param['sizeSubaperture'      ] = param['diameter']/param['nSubaperture']        # size of a sub-aperture projected in the M1 space
    param['samplingTime'         ] = 1/1000                                         # loop sampling time in [s]
    param['centralObstruction'   ] = 0.115                                              # central obstruction in percentage of the diameter
    param['nMissingSegments'     ] = 0                                             # number of missing segments on the M1 pupil
    param['m1_reflectivity'      ] = 1                                   # reflectivity of the 798 segments

    param['pixelSize'            ] = param['diameter']/param['resolution']          # number of missing segments on the M1 pupil

          
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NGS PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['magnitude'            ] = 8                                              # magnitude of the guide star
    param['opticalBand'          ] = 'I'                                            # optical band of the guide star
    
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param['nActuator'            ] = 676                                             # number of actuators 
    param['mechanicalCoupling'   ] = 0.45
    # mis-registrations                                                             
    param['shiftX'               ] = 0                                              # shift X of the DM in pixel size units ( tel.D/tel.resolution ) 
    param['shiftY'               ] = 0                                              # shift Y of the DM in pixel size units ( tel.D/tel.resolution )
    param['rotationAngle'        ] = 0                                              # rotation angle of the DM in [degrees]
    param['anamorphosisAngle'    ] = 0                                              # anamorphosis angle of the DM in [degrees]
    param['radialScaling'        ] = 0                                              # radial scaling in percentage of diameter
    param['tangentialScaling'    ] = 0                                              # tangential scaling in percentage of diameter
    param['isM4'                 ] = False                                             # number of actuators 
    param['isLBT']               = True
    param['pitch']               = 0.3                                              # pitch in metres    
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WFS PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['modulation'            ] = 3                                             # modulation radius in ratio of wavelength over telescope diameter
    param['pupilSeparationRatio'  ] = 1.2                                           # separation ratio between the PWFS pupils
    param['psfCentering'          ] = False                                          # centering of the FFT and of the PWFS mask on the 4 central pixels
    param['calibrationModulation' ] = 50                                            # modulation radius used to select the valid pixels
    param['lightThreshold'        ] = 0.1                                           # light threshold to select the valid pixels
    param['edgePixel'             ] = 1                                             # number of pixel on the external edge of the PWFS pupils
    param['extraModulationFactor' ] = 0                                             # factor to add/remove 4 modulation points (one for each PWFS face)
    param['postProcessing'        ] = 'slopesMaps'                                   # post-processing of the PWFS signals
    param['unitCalibration'       ] = False                                         # calibration of the PWFS units using a ramp of Tip/Tilt    
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP AND CALIBRATION DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param['nLoop'                 ] = 500                                           # number of iteration                             
    param['photonNoise'           ] = True
    param['readoutNoise'          ] = 0
    param['gainCL'                ] = 0.5
    param['nModes'                ] = 300
    param['nPhotonPerSubaperture' ] = 1000  
    param['getProjector'          ] = True
    
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # name of the system
    param['name'] = 'LBT_' +  param['opticalBand'] +'_band_'+ str(param['nSubaperture'])+'x'+ str(param['nSubaperture'])  
    # location of the calibration data
    param['pathInput'            ] = 'C:/Users/tsl29789/OneDrive - Science and Technology Facilities Council/Documents/LBT/SPRINT@LBT/psim/data_calibration/' 
    # location of the output data
    param['pathOutput'            ] = 'C:/Users/tsl29789/OneDrive - Science and Technology Facilities Council/Documents/LBT/SPRINT@LBT/psim/data_cl'
    

    print('Reading/Writting calibration data from ' + param['pathInput'])
    print('Writing output data in ' + param['pathOutput'])

    createFolder(param['pathOutput'])

    return param