#!/home/guillermo.valdes/.conda/envs/hht-py37/bin/python
# 
# Copyright (C) LIGO Scientific Collaboration (2019-)
#
# This code requests the data from a given channel, extracts the scattering 
# information using Hilbert-Huang transform, then estimates the fringe 
# frequency of scattered light based on the top-mass motion of LIGO suspensions,
# and computes the correlation between the scattering and the estimated fringes.
#
# Guillermo Valdes / Apr 2019
# 
# Usage:
# ./CheckScatCorr.py observatory gps -p primary-channel -d duration -f suspension-channel-list
# 
# Example 1:
# ./CheckScatCorr.py L1 1129835218
#
# Example 2:
# ./CheckScatCorr.py L1 1129835218 -p GDS-CALIB_STRAIN -d 100 -f Channel_List.txt
# output: Maximum correlation of 0.50 with channel L1:SUS-SRM_M1_DAMP_L_IN1_DQ

import os
import argparse

from gwpy.timeseries import TimeSeries
from gwpy.signal.filter_design import (notch, highpass, lowpass, concatenate_zpks)

from scipy.signal import tukey
from scipy.signal import butter, lfilter, freqz

import pylab as plt
import matplotlib as mpl; mpl.use('Agg')
import numpy  as np

from PyEMD import EMD
from scipy.signal import hilbert

# command line options
parser = argparse.ArgumentParser()
parser.add_argument("ifo", help="ifo name")
parser.add_argument("gps", help="centered gps time", type=float)
parser.add_argument('-p', '--channel', help="channel name", default='GDS-CALIB_STRAIN')
parser.add_argument('-d', '--duration', help="duration of segment", type=float, default=100)
parser.add_argument('-f', '--channel-file', help='path for channel file', type=os.path.abspath, default='Channel_List.txt')
args = parser.parse_args()

print(args)

# define parameters
ifo = args.ifo
chname = '%s:%s' % (ifo,args.channel)
duration = args.duration
gps  = args.gps
gpse = gps + duration/2
gpss = gps - duration/2

# load timeseries for witness channel
data = TimeSeries.get(chname, gpss, gpse)

# notch filters
n1 = notch(7.6, 1/data.dt)
n2 = notch(9.8, 1/data.dt)
n3 = notch(11, 1/data.dt)
n4 = notch(13.6, 1/data.dt)
n5 = notch(25, 1/data.dt)
n6 = notch(36.0, 1/data.dt)

zpk = concatenate_zpks(n1, n2, n3, n4, n5, n6)
filtered = data.filter(zpk, filtfilt=True)

# pre-processing witness channel data
bp = filtered.bandpass(15, 30)          # band-pass filter
white = bp.whiten(highpass=15)          # whitining filter
cdata = white[163840:1474560]           # crop data
window = tukey(len(cdata),alpha=0.25)   # tukey window
wdata = cdata*window                    # force extremes to zero
wdata = wdata.resample(256)             # decimate to 256 Hz sampling frequency

# time and data arrays
time_array = wdata.times - wdata.times[0]
signal = wdata.value

# empirical mode decomposition
print('Applying emd')
emd = EMD()
imf = emd.emd(signal)
print('Finish emd')

# hilbert spectral analysis
analytic_signal = hilbert(imf)
amplitude = np.abs(analytic_signal)
phase = np.unwrap(np.angle(analytic_signal))
instant_freq = np.diff(phase)/(2*np.pi*wdata.dt)

# processing instantaneous amplitude
order = 2               # filter order
fs = 1/wdata.dt         # sample rate, Hz
cutoff = 5.0            # desired cutoff frequency of the filter, Hz

b, a = butter(order, cutoff, btype='low', fs=256)
filtered_ia = lfilter(b, a, amplitude[0])

# call suspension channels
with open(args.channel_file, 'r') as f:
    channels = [name.rstrip('\n') for name in f]

# scattering predictor function     
def ScatteringPredictor(chname):
    position = TimeSeries.get(chname, gpss, gpse, verbose=True)
    position = position.bandpass(0.03, 10) # lowpass position channel
    velocity = (position[1:]-position[:-1])*(position.sample_rate.value)
    fudge = 1.0
    scatf1 = abs(fudge*8.0*velocity.value/1.064) 
    times = velocity.times - velocity.times[0]
    return scatf1, times

# estimate scattering predictors and correlation
l = []
p = []
for channel in channels:
    position_chan = '%s:%s' % (ifo,channel)
    scatf1, times = ScatteringPredictor(position_chan)
    print(position_chan)    

    # correlation estimation
    cc = np.corrcoef(scatf1[2559:-2560],filtered_ia)
    print(cc[0][1])
    corr_coef = ("%.2f" % cc[0][1])
    l.append(cc[0][1])
    p.append(position_chan)
    
# maximum correlation
mx_ind = np.argmax(l)    
max_corr = ("%.2f" % l[mx_ind])
max_chan = (p[mx_ind])

# print results
winner_corr = '\nMaximum correlation of %s ' % max_corr 
winner_chan = 'with channel %s\n' % max_chan
print(winner_corr + winner_chan)


