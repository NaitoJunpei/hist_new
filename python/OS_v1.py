##########
# OS_v1.py computes the optimal number of bins of time-histogram based on the optimization method proposed by Omi and Shinomoto, which may be applicable to non-Poisson spike trains. 
# needs libraries: (matplotlib, numpy, pandas). 

# Instruction
# put OS_v1.py in a folder. 
# import OS_v1
# then you may obtain OS_v1.().

# you need only OS function.
# the function OS takes a spike train as an argument. 
# spike train could be given by a list or numpy.array.
# the program selects the optimal bin size for a given spike train and draws the histogram. 
# references:
#  Takahiro Omi & Shigeru Shinomoto, "Optimizing time histograms for non-Poissonian spike trains", Neural Computation 23, 3125 (2011).
# Contact:
# Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp

##########

import matplotlib.pyplot as plt
import numpy as np

def OS(spike_times) :
    spike_times = np.array(spike_times)
    max_value   = max(spike_times)
    min_value   = min(spike_times)
    onset       = min_value - 0.001 * (max_value - min_value)
    offset      = max_value + 0.001 * (max_value - min_value)
    lv          = 0
    ISI         = np.diff(spike_times)

    for i in range(0, len(spike_times) - 2) :
        interval1 = ISI[i]
        interval2 = ISI[i + 1]

        if(interval1 + interval2 != 0) :
            lv += 3 * pow(interval1 - interval2, 2) / (pow(interval1 + interval2, 2) * (len(spike_times) - 2))
        else :
            lv += 3 / (len(spike_times) - 2)

    for bin_num in range(1, 500) :
        times = 10
        cost = cost_av(spike_times, onset, offset, lv, bin_num, times)
        
        if (bin_num == 1 or cost < cost_min) :
            cost_min        = cost
            optimal_bin_num = bin_num

    drawOS(spike_times, optimal_bin_num)

########## 
# cost_f 
# computes the cost function defined by Omi and Shinomoto

# arguments:
# spike_times: spike train
# start: time of the initial spike
# end: time of the final spike
# lv: the value of local variation Lv, which measures the spiking irregularity 
# bin_num: number of bins

# returns the cost function
########## 


def cost_f(spike_times, start, end, lv, bin_num) :
    bin_width = (end - start) / bin_num
    hist = np.histogram(spike_times, np.linspace(start, end, bin_num + 1))[0]

    fano = 2.0 * lv / (3.0 - lv)

    av   = np.mean(hist)
    va   = np.mean(hist * hist)
    w_av = np.mean(hist * fano)
    fano_bin = np.where(hist > 2, fano, 1.0)

    return ((2.0 * np.mean(hist * fano_bin) - (va - av * av)) / (bin_width * bin_width))

########## 
# cost_av
# computes an average cost function with respect to initial binning positions.

# arguments:
# spike_times: spike train
# onset: time of an initial spike
# offset: time of a final spike
# lv: the value of local variation Lv, which measures the spiking irregularity
# bin_num: the number of bins
# times: the number of initial binning positions

# returns the averaged cost function
##########


def cost_av(spike_times, onset, offset, lv, bin_num, times) :
    temp = 0.0
    bin_width = (offset - onset) / bin_num
    TT = np.hstack([spike_times, spike_times + (offset - onset)])
    for i in range(0, times) :
        start = onset + i * bin_width / times
        end = offset + i * bin_width / times
        temp += cost_f(TT, start, end, lv, bin_num)

    return temp / times

##########
# drawOS 
# draws a histogram

# arguments:
# spike_times: a spike train
# optimal_bin_num: an optimal number of bins
########## 

def drawOS(spike_times, optimal_bin_num):
    plt.hist(spike_times, optimal_bin_num)
    plt.yticks([])
    plt.show()
