# -*- coding: utf-8 -*-
""" 
CO TDLAS NO FITTING JUST FOR PLOTTING ON THE FLY 

Created on Tue Aug 11 14:36:13 2026

@author: ezb0082

aug11: establish file reading, signal isolating, etcetc
aug12: read asc file, pick latest sig file, fix shot count input, resoltion to conv to volts within sig

PICKING LATEST FILE ONLY WORKS FOR .SIG!!!!! 

"""
shotCount = 500 # CHANGE THIS TO HOW MANY SHOTS YOU TAKE !!!!

# CHANGE THIS TO THE FOLDER LOCATION
folderName = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\test folder" 

#%% packages
import os 
import numpy as np
import matplotlib.pyplot as plt

from pybaselines import Baseline
from pathlib import Path
from random import randint
#%% ----------------------------------------------------------------------------------------------------------------------
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

# signal = np.mean(signal_list, axis=0) # calc mean of all frames
#%% Read data as txt PRECONVERTED FROM GAGE .SIG FILE 
# start = time.time()
# folderName = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\Ramp\\CO Cell"

# signal_list = [] # allocate storage
# for frame in range(1,shotCount+1):
#     fName = "CO_42.75degC_height_5mm_1kHz_ramp_COcell_Record{frame:03d}.txt".format(frame=frame) # data file names
#     fPath = Path(folderName, fName) # combine path
#     rawData = np.genfromtxt(fPath, skip_header=13) # read data
#     signal_list.append(rawData) # append to d1ata list
# signal_list = np.array(signal_list) # append to array

# signal = np.mean(signal_list, axis=0) # calc mean of all 500 frames
# end = time.time()

# elapsed = end - start

# print(elapsed)

# plt.plot(signal)

#%% read data from .sig SPECIFIED

# f = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\Ramp\\CO_42.75degC_height_5mm_1kHz_ramp_phi_2.sig"
# dat = np.memmap(f, np.int16, 'r')
# header = dat[0:256]

# resolution = abs(int(header[150])) ## this is the factor to convert to volts !! 

# realData = dat[256:len(dat)] / resolution
# signal_list = np.reshape(realData, (shotCount, -1))

# signal = np.mean(signal_list, axis=0)
# signal = signal 

# plt.plot(signal)

#%% read data from .asc VOLTS 

# f = "C:\\Users\\ezb0082\\Downloads\\2026.08.11\\testGageScopeSaveVolts03.asc"

# dat = np.loadtxt(f, skiprows=13)
# signal_list = np.reshape(dat, (shotCount, -1))
# signal = np.mean(signal_list, axis=0)

# plt.plot(signal)

#%% ---------------------------------------------------------------------------------------------------------------------
#%% pick and plot 10 frames

fig, axs = plt.subplots(5, 2, sharex=True, figsize=[10,6]) # initialize figure
fig.supylabel("Voltage")  
fig.supxlabel("Nanoseconds")
plt.subplots_adjust(hspace=0.5) # space subplots apart

randShot = [] # allcoate stroage
for i in range(10): 
    randNumber = randint(1,500) # pick random frame
    randShot.append(randNumber) # append to array
randShot = np.array(randShot)     # make array

randShot.sort() # put random picks in order

for i, ax in zip(range(10), axs.ravel()):
    shotNumber = randShot[i] # pull shot number
    randSignal = signal_list[shotNumber,:] # pull shot data
    
    ax.plot(randSignal) # plot this signal
    ax.title.set_text('Shot {:}'.format(str(shotNumber))) # plot title

#%% isolate just the signal 

x = np.arange(0, len(signal), 1) # create time arrays for x axis 
x_mask = (x >= 8000) & (x <= 63900) # isloate the roi, the ramp function
x_narrow = x[x_mask]  # apply the mask and convert from ns to s 

sigMasked = signal[x_mask] # isolate the signal 

#%% remove the background

bg_fit_mask = ((x_narrow <= 31000) | (x_narrow >= 44500)) # implement mast  
# CO CELL FEATURE MASK: ((x_narrow <= 31000) | (x_narrow >= 44500)) 
# HENCKEN BURNER CO FEATURE: ((x_narrow <= 33750) | (x_narrow >= 44500)) 

baseline_fitter = Baseline(x_data=x_narrow) # init baseline fitter
fit, params_2 = baseline_fitter.asls(sigMasked, lam=1e6, p=0.001) # fit asls baseline to data
fitted_BG = np.interp(x_narrow, x_narrow[bg_fit_mask], fit[bg_fit_mask]) # apply mask to fit line and interp back to regular xaxis

#%% plot raw, isolated, and bg fitting 

fig, axs = plt.subplots(1, 3, sharey=True, figsize=[10,6]) # initialize figure
fig.supylabel("Voltage")
fig.supxlabel("Nanoseconds")

axs[0].plot(x, signal) # plot raw data
axs[0].title.set_text('Raw Data')

axs[1].plot(x_narrow, sigMasked) # plot masked data
axs[1].title.set_text('Isolated Signal')

axs[2].plot(x_narrow, sigMasked, label='Signal') # plot the isolated signal
axs[2].plot(x_narrow, fitted_BG, label="Background") # plot the fited backgroiund 
axs[2].title.set_text('Background Fitting')
axs[2].legend()


