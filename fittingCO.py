# -*- coding: utf-8 -*-
""" 
CO TDLAS FITTING FOR TEMPERATURE AND MOLE FRACTION OF CO 

Created on Tue Aug 11 14:36:13 2026

@author: ezb0082

aug12: establish fitting routine based on data reading code  
""" 

# CHANGE THIS TO HOW MANY SAMPLES YOU TAKE !!!!
shotCount = 500 

# CHANGE THIS TO THE FOLDER LOCATION
folderName = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\test folder" 

# This can be changed to better guess T and X
T = 2000 # k
X = 0.1 # CO mole fraction

# change this to your experimental path length
length = 2.54 # centimeters
passes = 1 # number of passes for multipass

#%% packages
import os 
import numpy as np
import matplotlib.pyplot as plt
import hapi as hp

from scipy.optimize import minimize, Bounds
from pybaselines import Baseline
from pathlib import Path
from scipy.io import loadmat
from scipy.signal import find_peaks

#%% functions

def hapi_calculation(species,p, t, x, length, resltn, afwing, wnb, dnu):
            nu, coef = hp.absorptionCoefficient_Voigt(SourceTables=species,    # calculate voigt absorption coefficient
                                              WavenumberRange=(wnb),
                                              WavenumberStep=(dnu),
                                              Environment={'p':p, 'T':t}, 
                                              Diluent={'Air':(1-x)},
                                              HITRAN_units=False)   
            coef_corrected = coef * x  # scale coefficent by mole fraction
            nu, absorp = hp.absorptionSpectrum(nu,coef_corrected,  # calculate absorption from abs coeff                 
                                              Environment={'l':length})  
            wl_convolve = 1e7/nu # calculate wavelength 
            wl_hapi, absorp_, i1, i2, slit = hp.convolveSpectrum(np.flip(wl_convolve),absorp, # convolve to insturment function   
                                              SlitFunction=hp.SLIT_GAUSSIAN,
                                              Resolution=resltn,AF_wing=afwing)
            wl = np.flip(wl_hapi) # wavelength
            
            return wl, absorp_ 

def objective(arguments, A_exp_wlc):
      
    # T, P, X, L = arguments[0], arguments[1], arguments[2], arguments[3] # pull t p x and l from the input array
    T = arguments[0] # pull t from the input array
    X = arguments[1] # pull t from the input array
    
    wl, A_sim = hapi_calculation('CO', P, T, X, length, resltn, afwing, wnb, dnu) # run Hcl calc for current condition

    A_sim_c = np.interp(wlc, np.flip(wl), np.flip(A_sim)) # interpolate to common wavelength scale`

    ssr = sum((A_exp_wlc - A_sim_c)** 2) # calculate r2

    return ssr

def wl_axis_cal(param, waveAxis, simAxis, abs_sim, abs_data):
    a1 = param[0]
 
    waveAxis_shift =  waveAxis + a1
    abs_sim_shift = np.interp(waveAxis_shift, simAxis, abs_sim)

    residuals = abs_data - abs_sim_shift
    residualrmse = np.sqrt(np.sum((residuals ** 2) / len(residuals)))  # root mean squared error
 
    return residualrmse
#%% ------------------------------------------------------------------------------------------------------------------
#%% pick latest SIG file in folder CHANGE FOLDERNAME TO location !!!
most_recent_file = None # initiallize
most_recent_time = 0  # initialize

for entry in os.scandir(folderName): # loop through folder
    if entry.is_file(): # check file
        mod_time = entry.stat().st_mtime_ns # ch3eck time
        if mod_time > most_recent_time: # compare tiems
            most_recent_file = entry.name # update file name
            most_recent_time = mod_time # update newest time
        
fullPath = Path(folderName, most_recent_file) # combine path
dat = np.memmap(fullPath, np.int16, 'r') # read data
header = dat[0:256] # isolate header

resolution = abs(int(header[150])) ## this is the factor to convert to volts !! 

realData = dat[256:len(dat)] / resolution # get real data 

signal_list = np.reshape(realData, (shotCount, -1)) # reshape array, unstack all shots from being vertically stacked

signal = np.mean(signal_list, axis=0) # take average of all shots

#%% Read data as txt PRECONVERTED FROM GAGE .SIG FILE 

# folderName = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\Ramp\\CO Cell"

# signal_list = [] # allocate storage
# for frame in range(1,shotCount+1):
#     fName = "CO_42.75degC_height_5mm_1kHz_ramp_COcell_Record{frame:03d}.txt".format(frame=frame) # data file names
#     fPath = Path(folderName, fName) # combine path
#     rawData = np.genfromtxt(fPath, skip_header=13) # read data
#     signal_list.append(rawData) # append to d1ata list
# signal_list = np.array(signal_list) # append to array

# signal = np.mean(signal_list, axis=0) # calc mean of all 500 frames


#%% read data from .sig

# f = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\Ramp\\CO_42.75degC_height_5mm_1kHz_ramp_phi_2.sig"
# dat = np.memmap(f, np.int16, 'r')
# header = dat[0:256]

# resolution = abs(int(header[150])) ## this is the factor to convert to volts !! 

# realData = dat[256:len(dat)] / resolution
# signal_list = np.reshape(realData, (shotCount, -1))

# signal = np.mean(signal_list, axis=0)

#%% read data from .asc VOLTS 

# f = "C:\\Users\\ezb0082\\Downloads\\2026.08.11\\testGageScopeSaveVolts03.asc"

# dat = np.loadtxt(f, skiprows=13)
# signal_list = np.reshape(dat, (shotCount, -1))
# signal = np.mean(signal_list, axis=0)

#%% Read data from .mat 
signal_list = []
folderName = "E:\\Elisabeth\\tdlas\\2026.8.14"
for frame in range(1,shotCount+1):
    fName = "2026.8.14_{frame:04d}.mat".format(frame=frame)
    fPath = Path(folderName, fName)
    dat = loadmat(fPath)
    datSignal = dat["A"]
    datReduced = datSignal[:,0]
    signal_list.append(datReduced)
signal_list = np.array(signal_list)
signal = np.mean(signal_list, axis=0)

bg_list = []
folderName = "E:\\Elisabeth\\tdlas\\2026.8.14-wo_flame"
for frame in range(1,500+1):
    fName = "2026.8.14-wo_flame_{frame:03d}.mat".format(frame=frame)
    fPath = Path(folderName, fName)
    bgdat = loadmat(fPath)
    bgdatSignal = bgdat["A"]
    bgdatReduced = bgdatSignal[:,0]
    bg_list.append(bgdatReduced)
bg_list = np.array(bg_list)
bg = np.mean(bg_list, axis=0)

#%% ---------------------------------------------------------------------------------------------------------------------
#%% isolate signal 
signal = signal - signal[0] # make starting pt ze4ro
bg = bg - bg[0] # make starting point zero 

signal_norm = signal.flatten() / np.max(signal) # normalize signal
bg_norm = bg.flatten() / np.max(bg) # normalize bg

x = np.arange(0, len(signal), 1) # create time arrays for x axis 

endRamp = find_peaks(signal, height= 0.1, distance = 50)
endRamp = endRamp[0]
x_mask = (x >= 2570) & (x <= endRamp[-1]) # isloate the roi
x_narrow = x[x_mask]  # apply the mask and convert from ns to s 

sigMasked = signal_norm[x_mask] # isolate the signal 
bg_masked = bg_norm[x_mask]

#%% transmission and absorption 

transmission = sigMasked / bg_masked # transmission I / Io 

absorption = 1 - transmission # absortion 1- transmission

fig, axs = plt.subplots(1, 2, figsize=[10,6])
fig.supxlabel("Time (ns)")

axs[0].plot(transmission)
axs[0].title.set_text('Transmission')

axs[1].plot(absorption)
axs[1].title.set_text('Absorption')

plt.show()

#%% denoise - CO ONY 
def denoise(x, y):
    baseline_fitter = Baseline(x_data=x)
    line, params_2 = baseline_fitter.drpls(y, lam=1e6)
    return line

absorption_DN = denoise(x_narrow, absorption)

# plt.figure()
# plt.plot(x_narrow, absorption)
# plt.plot(x_narrow, absorption_DN)

absorption = absorption_DN
#%% baseline correct 

deg = 2
poly = np.polyfit(x_narrow, absorption, deg)
bkg = np.polyval(poly, x_narrow)

absorption_C = absorption - bkg 

# plt.figure()
# plt.plot(x_narrow, absorption)
# plt.plot(x_narrow, bkg)
# plt.plot(x_narrow, absorption_C)
# plt.show()

#%% temperature and mole fraction fitting

sort_idx = np.argsort(x_narrow) # sort the x axis 
xs = x_narrow[sort_idx] 

wavelength_exp = 2329.7 + 0.95 * (xs - xs.min()) / (xs.max() - xs.min()) # build wl axis initial gues

hp.db_begin('HAPI_DATA') # initialize hapi database

total_pathLength = length * passes  # cm
resltn = 0.0001 # cm-1
afwing = resltn * 5  # instrument function width

wnb = [4290, 4293] # ~ 2329 nm to 2330 nm CO feature range
dnu = 0.00001 # wavenumber stepsize, must be less than resolution

P = 1 # PRessure in atm 

wlc = np.arange(2329.95, 2330.5, 0.001) # establish common wavelengthj sacle 

wl, absorp_sim = hapi_calculation('CO', P, T, X, length, resltn, afwing, wnb, dnu) # calculate absorp. for CEA values
    
A_sim_wlc = np.interp(wlc, np.flip(wl), np.flip(absorp_sim)) # interpolate simulated spectrra to common scale

bnds = Bounds([-1e6], [1e6]) # bounds for wl axis fitting
init = [1e-3] # initial guess for wl axis fitting
wl_fit_result = minimize(wl_axis_cal, init, bounds=bnds, tol= 1e-6, method='Nelder-Mead', args=(wavelength_exp, wlc, A_sim_wlc, absorption_C)) # fit the wavelength axis 
wl_axis_fit =  wavelength_exp + wl_fit_result.x[0] # adjust wl axis

A_exp_wlc = np.interp(wlc, wl_axis_fit, absorption_C) # interpolate experimental spectra to common scale using fitted wl axis

bnds = Bounds([200, 0.01], # lower bounds
              [4000, 1]) # upper bounds
    
initial_guess = [T, X ] # initial guess for flame  condtions

result = minimize(objective, initial_guess, bounds=bnds, tol=1e-4, method='Nelder-Mead', args=(A_exp_wlc)) # run minimizing routine for  speactrra fittig

print(result) # print minimizing routine results

T_fit =  result.x[0] # read result for T
X_fit = result.x[1] # read result for X

wl_min, A_minimized = hapi_calculation('CO', P, T_fit, X_fit, length, resltn, afwing, wnb, dnu)  # calculate absorb at fitted t and x

A_min_wlc = np.interp(wlc, np.flip(wl_min), np.flip(A_minimized)) # interpolate minimized to common scale
    
plt.figure() # plot
plt.plot(wlc, A_exp_wlc , label='Exp.')
plt.plot(wlc, A_min_wlc, label='Fit T = {:.2f} K, X = {:.3f}'.format(T_fit, X_fit))
plt.legend()
plt.title("Experimental and Fitted Spectra")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorption")
plt.show()
