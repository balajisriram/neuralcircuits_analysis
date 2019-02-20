# convert_pl2_to_dat.py - PyPL2 usage example
#
# (c) 2016 Plexon, Inc., Dallas, Texas
# www.plexon.com - support@plexon.com
#
# This software is provided as-is, without any warranty.
# You are free to modify or share this file, provided that the above
# copyright notice is kept intact.
from __future__ import print_function
from pypl2 import pl2_ad, pl2_spikes, pl2_events, pl2_info
import sys
import os
import numpy as np
import scipy.signal as sig
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def load_and_pack(filename):
    spkinfo, evtinfo, adinfo = pl2_info(filename)
    print("\nContinuous A/D Channel Info from pl2_info()")
    print("\n# Channel Name\tCount")
    print("- ------------\t-----")
    wideband_signals = []
    n_samples = []
    for n in range(len(adinfo)):
        if 'WB' in adinfo[n].name and adinfo[n].n>0:
            print("%s %-10s\t%s" % (adinfo[n].channel, adinfo[n].name, adinfo[n].n))
            wideband_signals.append(adinfo[n].name)
            if not n_samples:
                n_samples = adinfo[n].n
            elif n_samples != adinfo[n].n:
                ValueError('a channel does not have the same number of data points as another')
            else:
                pass #all ok
    print("- ------------\t-----")
    n_channels = len(wideband_signals)
    print("Found %d WB signals with %d points each. Going to numpy and pack" % (n_channels, n_samples))
    data_array = np.zeros([n_samples, n_channels], np.int16)

    print('finding max absolute value across all channels')
    max_value = np.NINF
    for n in range(n_channels):
        print('%d.' %n ,end='')
        sys.stdout.flush()
        temp1 = pl2_ad(filename, n)
        temp = np.asarray(temp1.ad)
        if np.max(np.abs(temp))>max_value:
            max_value = np.max(np.abs(temp))
    print("")
    print("maximum absolute value found = %s" % max_value)

    print("loading, normalizing and packing")
    for n in range(n_channels):
        print('%d.' %n ,end='')
        sys.stdout.flush()
        temp1 = pl2_ad(filename, n)
        temp = np.asarray(temp1.ad)
        temp = np.int16(temp/max_value*32767)
        data_array[:,n] = temp
    print("")
    temp,junk = os.path.splitext(filename)
    dat_filename = temp+'.dat'
    print("Writing to %s" % dat_filename)
    data_array.tofile(dat_filename)

def get_running(filename):
    temp = pl2_ad(filename,'AIF01')
    running = np.asarray(temp.ad)
    running_rate = temp.adfrequency
    
    duration_for_median_filter = 0.001 # secs
    samples_for_med_filt = np.int(duration_for_median_filter*running_rate)
    
    filt_running = sig.medfilt(running, samples_for_med_filt)
    
    #plt.plot(running[0:10000],'k')
    # plt.plot(filt_running,'r')
    # plt.show()
    
    plt.plot(filt_running,'k')
    pctile_80 = np.percentile(filt_running,80)
    print('80th_percentile::',pctile_80)
    higher_than_pctile = filt_running
    higher_than_pctile[higher_than_pctile<=pctile_80] = np.nan
    plt.plot(higher_than_pctile,'g')
    plt.show()
    
def fourier1(x,a,b):
    return a*np.cos(2*np.pi*60.*x+b)
    
def get_info(filename):
    spkinfo, evtinfo, adinfo = pl2_info(filename)
    print("\nContinuous A/D Channel Info from pl2_info()")
    print("\n# Channel Name\tCount")
    print("- ------------\t-----")
    wideband_signals = []
    n_samples = []
    for n in range(len(adinfo)):
        # if 'WB' in adinfo[n].name:
        print("%s %-10s\t%s" % (adinfo[n].channel, adinfo[n].name, adinfo[n].n))
            # wideband_signals.append(adinfo[n].name)
            # if not n_samples:
                # n_samples = adinfo[n].n
            # elif n_samples != adinfo[n].n:
                # ValueError('a channel does not have the same number of data points as another')
            # else:
                # pass #all ok
                
    #Get event data on select channels and print out interesting information
    print("\nEvent Data from pl2_events()")
    print("\nEvent   Number of Events First Timestamp (s)")
    print("------- ---------------- -------------------")
    
    for n in range(len(evtinfo)):
        evt = pl2_events(filename, evtinfo[n].name)
        print("%-7s %-16s %s" % (evtinfo[n].name, evt.n, evt.timestamps[0]))

    temp1 = pl2_ad(filename, 0)
    print("sampling frequency = %s" % temp1.adfrequency)
    
def compare_WB_and_FP(filename):
    # get the data
    print('getting data')
    temp1 = pl2_ad(filename,0)
    wb1 = np.asarray(temp1.ad)
    wb1 = wb1[0:1200000]
    samp_freq = temp1.adfrequency
    wb2,filtered_60,fourier_fit = remove_line_noise(wb1,fs=samp_freq)
    
    # filter
    print('filtering')
    (b,a) = sig.butter(2,200./samp_freq, 'lowpass')
    filt1 = sig.filtfilt(b,a,wb2,method='pad')
    print('done')
    
    # plot
    plt.subplot(1,2,1)
    plt.plot(np.arange(1200000)/samp_freq,filtered_60,color='gray')
    # plt.plot(np.arange(1200000)/samp_freq,filtered_60-fourier_fit,color='black')
    plt.subplot(1,2,2)
    f1,density1=sig.periodogram(filt1,samp_freq)
    f2,density2=sig.periodogram(wb1,samp_freq)
    which1 = np.bitwise_and(f1>1,f1<200)
    which2 = np.bitwise_and(f2>1,f2<200)
    #plt.semilogy(f1[which1],density1[which1],color='black')
    plt.semilogy(f2[which2],density2[which2],color='gray')
    
    plt.show()
  
def remove_line_noise_bad(inp, fs=1, n_harmonics=1, line_f_base=60., err_tol=0.001, hw = 2, win_type='hamming'):
    # adapted from https://www.mathworks.com/matlabcentral/fileexchange/54228-remove-line-noise
    inp = np.double(inp)
    line_f = np.arange(1,n_harmonics+1)*line_f_base
    # get the window size
    z,err,id = 4,np.INF,-1
    while err/line_f[0] > err_tol or id<(hw*2-1):
        z += 1
        W = fs*np.linspace(0.,1.,2**z)
        err = np.min(np.abs(W-line_f[0]))
        id = np.argmin(np.abs(W-line_f[0]))
    m = 2**z
    
    W = np.linspace(0.,1.,m)
    line_id = np.zeros_like(line_f)
    
    for i,hz in enumerate(line_f):
        line_id[i] = np.argmin(np.abs(W-hz))
    
    n_step = 2
    step = np.round(m/n_step)
    
    print('Removing line noise by spectral estimation')
    
    if win_type=='hamming':
        window = sig.hamming(m,sym=False)
    elif win_type=='hanning':
        window = sig.hanning(m,sym=False)
        
    x = np.arange(m/2+1,m)
    y = inp[0:m/2]

def remove_line_noise(wb, fs=40000.):
    b,a = sig.butter(3,[55*2/fs, 65*2/fs],'bandpass',analog=True)
    filt = sig.filtfilt(b,a,wb,method='pad')
    t_max = wb.size/fs
    t = np.linspace(0,t_max,wb.size)
    popt, pcov = curve_fit(fourier1, t, filt)
    
    # plt.plot(t,wb,'black')
    plt.plot(t,filt,'gray')
    plt.plot(t,filt-fourier1(t,*popt),'red',alpha=0.5)
    plt.show()
    
    return wb-fourier1(t,*popt),filt,fourier1(t,*popt)
    
    
def _test_fourier_fit():
    t = np.linspace(0,10,400000)
    y = np.sin(2*np.pi*2*t+np.pi/4) + 0.3*np.sin(2*np.pi*8*t+np.pi/4) + 0.1*np.random.randn(400000)
    popt, pcov = curve_fit(fourier1, t, y)
    fit_val = fourier1(t,*popt)
    plt.plot(t,y)
    plt.plot(t,fit_val,'r')
    plt.plot(t,y-fit_val,'gray')
    
    plt.show()
    
        
if __name__ == "__main__":
    # load_and_pack(r'C:\Users\bsriram\Desktop\BiogenData\102618_1\File1.pl2')
    # get_info(r'C:\Users\bsriram\Desktop\BiogenData\102618_1\File1.pl2')
    get_running(r'C:\Users\bsriram\Desktop\BiogenData\102618_1\File1.pl2')
    # compare_WB_and_FP(r'C:\Users\bsriram\Desktop\BiogenData\102618_1\File1.pl2')
    # remove_line_noise()


    # #Get continuous a/d data on first channel and print out interesting information
    # ad = pl2_ad(filename, 0)
    # print "\nContinuous A/D Channel 0 Data from pl2_ad()"
    # print "\nFrequency Number of Points First Four A/D Points (mV)"
    # print "--------- ---------------- ---------------------"
    # print "%-10s%-17s%s, %s, %s, %s" % (int(ad.adfrequency),
    #                                     ad.n,
    #                                     ad.ad[0] * 1000,
    #                                     ad.ad[1] * 1000,
    #                                     ad.ad[2] * 1000,
    #                                     ad.ad[3] * 1000)
    #
    # #Get spikes on first channel and print out interesting information on the first four spikes.
    # spikes = pl2_spikes(filename, 0)
    # print "\nSpike Channel 0 Data for First Four Waveforms from pl2_spikes()"
    # print "\nTimestamps (s) Unit First Four Waveform Points (uV)"
    # print "-------------- ---- -------------------------------"
    # for n in range(4):
    #     print "%-15s%-5s%s, %s, %s, %s" % (spikes.timestamps[n],
    #                                        spikes.units[n],
    #                                        spikes.waveforms[n][0] * 1000000,
    #                                        spikes.waveforms[n][1] * 1000000,
    #                                        spikes.waveforms[n][2] * 1000000,
    #                                        spikes.waveforms[n][3] * 1000000)
    #
    # #Get event data on select channels and print out interesting information
    # print "\nEvent Data from pl2_events()"
    # print "\nEvent   Number of Events First Timestamp (s)"
    # print "------- ---------------- -------------------"
    #
    # for n in range(len(evtinfo)):
    #     evt = pl2_events(filename, evtinfo[n].name)
    #     print "%-7s %-16s %s" % (evtinfo[n].name, evt.n, evt.timestamps[0])
    #
    # #Get strobed event data and print out interesting information
    # strobedevt = pl2_events(filename, 'Strobed')
    # print "\nFirst Ten Strobe Values from pl2_events()"
    # print "\nStrobe Value Timestamp (s)"
    # print "------------ -------------"
    # for n in range(10):
    #     print ("%-12s %s") % (strobedevt.values[n], strobedevt.timestamps[n])
