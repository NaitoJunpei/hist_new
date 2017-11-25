##########
# SS_v2.py computes the optimal number of bins of time-histogram based on the optimization method proposed by Shimazaki and Shinomoto 2007. 
# needs libraries: (matplotlib, numpy).

# Instruction
# put SS_v2.py in a folder.
# import SS_v2
# then you may obtain SS_v2.().
 
# you need only SS function.
# the function SS takes a spike train as an argument.
# spike train could be given by list or numpy.array.
# the program selects the optimal bin size for a given spike train and draws the histogram. 
# references:
#  H. Shimazaki and S. Shinomoto, A method for selecting the bin size of a time histogram. Neural Computation (2007) 19:1503-1700. 
# Shigeru Shinomoto (2010) Estimating the firing rate. in "Analysis of Parallel Spike Train Data" (eds. S. Gruen and S. Rotter) (Springer, New York).
# Contact:
# Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp
########## 
# SS_v2.py
# Junpei Naito 2017/11/14
########## 

import matplotlib.pyplot as plt
import numpy as np

def SS(spike_times) :
    spike_times = np.array(spike_times)
    max_value   = max(spike_times)
    min_value   = min(spike_times)
    onset       = min_value - 0.001 * (max_value - min_value)
    offset      = max_value + 0.001 * (max_value - min_value)

    # computes the cost function by changing the number of bins
    # adopts the number of bins that minimizes the cost function
    
    for bin_num in range(1, 500) :
        cost = cost_av(spike_times, onset, offset, bin_num, 10)
        if (bin_num == 1 or cost < cost_min):
            cost_min        = cost
            optimal_bin_num = bin_num

    drawSS(spike_times, optimal_bin_num)
    return optimal_bin_num

########## 
# cost_f 
# computes the cost function defined by Shimazaki and Shinomoto

# arguments:
# spike_times: spike train
# start: time of the initial spike
# end: time of the final spike
# bin_num: number of bins

# returns the cost function 
########## 

def cost_f(spike_times, start, end, bin_num) :
    bin_width = (end - start) / bin_num
    hist = np.histogram(spike_times, np.linspace(start, end, bin_num + 1))[0]

    av   = np.mean(hist)
    va   = np.mean(hist * hist)

    return ((2.0 * av - (va - av * av)) / (bin_width * bin_width))

########## 
# cost_av
# computes an average cost function with respect to initial binning positions.

# arguments:
# spike_times: spike train
# onset: time of an initial spike
# offset: time of a final spike
# bin_num: the number of bins
# times: the number of initial binning positions

# returns the averaged cost function
##########

def cost_av(spike_times, onset, offset, bin_num, times) :
    temp = 0.0
    bin_width = (offset - onset) / bin_num
    TT = np.hstack([spike_times, spike_times + (offset - onset)])

    # averages the cost with respect to the starting positions.
    # times: number of starting positions.

    for i in range(0, times) :
        start = onset + i * bin_width / times
        end = offset + i * bin_width / times
        temp += cost_f(TT, start, end, bin_num)

    return temp / times

##########
# drawSS
# draws a histogram

# arguments:
# spike_times: a spike train
# optimal_bin_num: an optimal number of bins
########## 

def drawSS(spike_times, optimal_bin_num):
    plt.hist(spike_times, optimal_bin_num)
    plt.yticks([])
    plt.show()
