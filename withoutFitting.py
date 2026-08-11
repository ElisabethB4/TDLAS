# -*- coding: utf-8 -*-
""" 
CO TDLAS NO FITTING JUST FOR PLOTTING ON THE FLY 

Created on Tue Aug 11 14:36:13 2026

@author: ezb0082

aug11: establish file reading, signal isolating, etcetc 

TO DO: 
    -- pick up latest file in folder
    -- wavelength axis  
"""
import numpy as np
import matplotlib.pyplot as plt

from pybaselines import Baseline
from pathlib import Path
from random import randint

#%% pick latest file 
""" WORKING ON THSI STILL ....  """

# folderName = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\Ramp"

# most_recent_file = None
# most_recent_time = 0 

# for entry in os.scandir(folderName):
#     if entry.is_file():
#         mod_time = entry.stat().st_mtime_ns
#         if mod_time > most_recent_time:
#             most_recent_file = entry.name
#             most_recent_time = mod_time
            
# with open(most_recent_file, 'r') as f:
#     dat = np.memmap(f, np.int16, 'r')
    
#%% Read data as txt PRECONVERTED FROM GAGE .SIG FILE 

# folderName = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.29\\without\\CO_42.75degC_height_5mm_1kHz_phi_2" # cell data path

# signal_list = [] # allocate storage
# for frame in range(1,501):
#     fName = "phi_2_Record{frame:03d}.txt".format(frame=frame) # data file names
#     fPath = Path(folderName, fName) # combine path
#     rawData = np.genfromtxt(fPath, skip_header=13) # read data
#     signal_list.append(rawData) # append to d1ata list
# signal_list = np.array(signal_list) # append to array

# signal = np.mean(signal_list, axis=0) # calc mean of all 500 frames

#%% read data from .sig

f = "C:\\Users\\ezb0082\\OneDrive - Auburn University\\.RESEARCH\\CO TDLAS\\2026.07.16\\CO Hencken Burner\\Ramp\\CO_42.75degC_height_5mm_1kHz_ramp_COcell.sig"

dat = np.memmap(f, np.int16, 'r')
header = dat[0:256]
realData = dat[256:len(dat)]
lenSig = round(len(realData) / 500)
signal_list = np.reshape(realData, (500, -1))

signal = np.mean(signal_list, axis=0)

#%% pick and plot 10 frames

fig, axs = plt.subplots(5, 2, sharex=True, figsize=[10,6]) # initialize figure
fig.supylabel("Voltage") 
fig.supxlabel("Nanoseconds")
plt.subplots_adjust(hspace=0.5) # space subplots apart

randShot = []
for i in range(10):
    randNumber = randint(1,500) # pick random frame
    randShot.append(randNumber)
randShot = np.array(randShot)    

randShot.sort() # put random picks in order

for i, ax in zip(range(10), axs.ravel()):
    shotNumber = randShot[i]
    randSignal = signal_list[shotNumber,:] # pull shot data
    
    ax.plot(randSignal)
    ax.title.set_text('Shot {:}'.format(str(shotNumber)))

#%% isolate just the signal 

x = np.arange(0, len(signal), 1) # create time arrays for x axis 
x_mask = (x >= 8000) & (x <= 63900) # isloate the roi, the ramp function
x_narrow = x[x_mask]  # apply the mask and convert from ns to s 

sigMasked = signal[x_mask] # isolate the signal 

#%% remove the background

bg_fit_mask = ((x_narrow <= 31000) | (x_narrow >= 44500)) # creat a mask for the feature THIS IS CELL
# HENCKEN BURNER CO FEATURE: ((x_narrow <= 33750) | (x_narrow >= 44500)) # creat a mask for the feature

baseline_fitter = Baseline(x_data=x_narrow) # init baseline fitter
fit, params_2 = baseline_fitter.asls(sigMasked, lam=1e6, p=0.001) # fit asls baseline to data
fitted_BG = np.interp(x_narrow, x_narrow[bg_fit_mask], fit[bg_fit_mask]) # apply mask to fit line and interp back to regular xaxis

#%% plot raw, isolated, and bg fitting 

fig, axs = plt.subplots(1, 3, sharey=True, figsize=[10,6])
fig.supylabel("Voltage")
fig.supxlabel("Nanoseconds")

axs[0].plot(x, signal)
axs[0].title.set_text('Raw Data')

axs[1].plot(x_narrow, sigMasked)
axs[1].title.set_text('Isolated Signal')

axs[2].plot(x_narrow, sigMasked, label='Signal')
axs[2].plot(x_narrow, fitted_BG, label="Background")
axs[2].title.set_text('Background Fitting')
axs[2].legend()

#%% WAvelength Axis 

"""
this is a probelm for later
"""

#%% calculate absorption using the fitted background

transmission = sigMasked / fitted_BG # transmission I / Io 

absorption = 1 - transmission # absortion 1- transmission

fig, axs = plt.subplots(1, 2, sharey=True, figsize=[10,6])
axs[0].plot(transmission)
axs[0].title.set_text('Transmission')

axs[1].plot(absorption)
axs[1].title.set_text('Absorption')
