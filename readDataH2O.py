# -*- coding: utf-8 -*-
""" 
H2O TDLAS NO FITTING JUST FOR PLOTTING ON THE FLY 

Created on Tue Aug 11 14:36:13 2026

@author: ezb0082

PICKING LATEST FILE ONLY WORKS FOR .SIG!!!!! 

"""
shotCount = 1300 # CHANGE THIS TO HOW MANY SHOTS YOU TAKE !!!!

# CHANGE THIS TO THE FOLDER LOCATION
folderName = "E:\\Elisabeth\\tdlas\\2026.8.14"

#%% packages
import os 
import numpy as np
import matplotlib.pyplot as plt

from pybaselines import Baseline
from pathlib import Path
from random import randint
from scipy.io import loadmat
from scipy.signal import find_peaks
#%% ----------------------------------------------------------------------------------------------------------------------
#%% pick latest SIG file in folder CHANGE FOLDERNAME TO location !!!
# most_recent_file = None # initiallize
# most_recent_time = 0  # initialize

# for entry in os.scandir(folderName): # loop through folder
#     if entry.is_file(): # check file
#         mod_time = entry.stat().st_mtime_ns # ch3eck time
#         if mod_time > most_recent_time: # compare tiems
#             most_recent_file = entry.name # update file name
#             most_recent_time = mod_time # update newest time
        
# fullPath = Path(folderName, most_recent_file) # combine path
# dat = np.memmap(fullPath, np.int16, 'r') # read data
# header = dat[0:256] # isolate header

# resolution = abs(int(header[150])) ## this is the factor to convert to volts !! 

# realData = dat[256:len(dat)] / resolution # get real data 

# signal_list = np.reshape(realData, (shotCount, -1)) # reshape array, unstack all shots from being vertically stacked

# signal = np.mean(signal_list, axis=0) # take average of all shots

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

#%% Read data from .mat 
signal_list = []
folderName = "C:\\Users\\elisa\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.08.14\\2026.8.14"
for frame in range(1,shotCount+1):
    fName = "2026.8.14_{frame:04d}.mat".format(frame=frame)
    fPath = Path(folderName, fName)
    dat = loadmat(fPath)
    datSignal = dat["D"]
    datReduced = datSignal[:,0]
    signal_list.append(datReduced)
signal_list = np.array(signal_list)
signal = np.mean(signal_list, axis=0)

bg_list = []
folderName = "C:\\Users\\elisa\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.08.14\\2026.8.14-wo_flame"
for frame in range(1,500+1):
    fName = "2026.8.14-wo_flame_{frame:03d}.mat".format(frame=frame)
    fPath = Path(folderName, fName)
    bgdat = loadmat(fPath)
    bgdatSignal = bgdat["D"]
    bgdatReduced = bgdatSignal[:,0]
    bg_list.append(bgdatReduced)
bg_list = np.array(bg_list)
bg = np.mean(bg_list, axis=0)

#%% ---------------------------------------------------------------------------------------------------------------------
#%% pick and plot 10 frames

fig, axs = plt.subplots(5, 2, sharex=True, figsize=[10,6]) # initialize figure
fig.supylabel("Voltage")  
fig.supxlabel("Time (ns)")
plt.subplots_adjust(hspace=0.5) # space subplots apart
fig.suptitle("10 H2O Samples")

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

plt.figure()
plt.ylabel("Voltage")
plt.xlabel("Time (ns)")
plt.title("10 H2O Samples Overlaid")
for i in range(len(randShot)):
    plt.plot(signal_list[randShot[i],:])
#%% isolate signal 
signal = signal - signal[0] # make starting pt ze4ro
bg = bg - bg[0] # make starting point zero 

signal_norm = signal / np.sum(signal) # normalize signal
bg_norm =bg / np.sum(bg) # normalize bg

x = np.arange(0, len(signal), 1) # create time arrays for x axis 

endRamp = find_peaks(signal, height= 0.1, distance = 100)
endRamp = endRamp[0]
x_mask = (x >= 1000) & (x <= endRamp[-1]) # isloate the roi
x_narrow = x[x_mask]  # apply the mask and convert from ns to s 

sigMasked = signal_norm[x_mask] # isolate the signal 
bg_masked = bg_norm[x_mask]
# sigMasked = sigMasked / np.mean(sigMasked)

#%% plot raw & isolated signal

fig, axs = plt.subplots(1, 2, sharey=True, figsize=[10,6]) # initialize figure
fig.supylabel("Normalized Signal")
fig.supxlabel("Time (ns)")
fig.suptitle("H2O, Averaged {:.0f} Shots".format(shotCount))

axs[0].plot(signal_norm) # plot raw data
axs[0].plot(bg_norm)
axs[0].title.set_text('Raw Data')

axs[1].plot(x_narrow, sigMasked) # plot masked data
axs[1].plot(x_narrow, bg_masked)
axs[1].title.set_text('Isolated Signal')

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
