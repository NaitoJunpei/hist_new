##########
# BRE.py returns Bayesian estimation of the rate of event occurrence.
# needs libraries: (matplotlib, numpy).

# Instruction
# put BRE.py on a folder in a path.
# import BRE
# then you may obtain BRE.().
# you need only BRE function
# the function BRE(Bayesian Rate Estimate) takes a spike train as an argument.
# spike train could be given by list or numpy.array.
# a parameter beta is determined by EM algorithm, and a figure is drawn from the rate estimated with the Kalman filter. 
# references:
#  S. Koyama and S. Shinomoto, Empirical Bayes interpretations of random point events.  J. Phys. A (2005) 38:L531-L537.
# T. Shimokawa and S. Shinomoto, Estimating instantaneous irregularity of neuronal firing. Neural Computation (2009) 21:1931-1951.
# Shigeru Shinomoto (2010) Estimating the firing rate. in "Analysis of Parallel Spike Train Data" (eds. S. Gruen and S. Rotter) (Springer, New York).
# Contact:
# Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp
##########

import matplotlib.pyplot as plt
import numpy as np
import math

def BRE(spike_times) :
    spike_times = np.array(list(spike_times))
    max_value   = max(spike_times)
    min_value   = min(spike_times)

    ISI    = np.diff(spike_times)
    mu     = len(spike_times) / (max_value - min_value)
    beta0  = pow(mu, -3)
    beta   = EMmethod(ISI, beta0)
    kalman = KalmanFilter(ISI, beta)

    drawBRE(spike_times, kalman)

####
# EMmethod
# estimates a parameter beta with the EM algorithm.

# arguments:
# ISI: inter spike interval
# beta0: initial value of the parameter beta

# returns the estimated beta.
####

def EMmethod(ISI, beta0) :
    N = len(ISI)
    beta = 0
    beta_new = beta0

    for j in range(0, 100) :
        beta = beta_new
        kalman = KalmanFilter(ISI, beta)

        beta_new = 0
        t0 = 0

        for i in range(0, N - 1) :
            if(ISI[i] > 0) :
                beta_new += (kalman[1][i + 1] + kalman[1][i] - 2 * kalman[2][i]
                             + (kalman[0][i + 1] - kalman[0][i])
                             * (kalman[0][i + 1] - kalman[0][i])) / ISI[i]
            else :
                t0 += 1

        beta_new = (N - t0 - 1) / (2 * beta_new)

    return beta_new

####
# KalmanFilter
# estimates the rate of event occurrence with the Kalman filtering.
# arguments:
# ISI: inter spike interval
# beta: the parameter
# returns the rate of event occurrence.
####

def KalmanFilter(ISI, beta) :
    N = len(ISI)
    IEL = N / sum(ISI)
    IVL = pow(IEL / 3, 2)
    A = IEL - ISI[0] * IVL
    EL = np.empty([2, N])
    VL = np.empty([2, N])

    EL_N = np.empty(N)
    VL_N = np.empty(N)
    COVL_N = np.empty(N)

    EL[0][0] = (A + math.sqrt(A * A + 4 * IVL)) / 2
    VL[0][0] = 1 / (1 / IVL + 1 / pow(EL[0][0], 2))

    # prediction and filtering
    for i in range(0, N - 1) :
        EL[1][i] = EL[0][i]
        VL[1][i] = VL[0][i] + ISI[i] / (2 * beta)

        A = EL[1][i] - ISI[i + 1] * VL[1][i]
        EL[0][i + 1] = (A + math.sqrt(A * A + 4 * VL[1][i])) / 2
        VL[0][i + 1] = 1 / (1 / VL[1][i] + 1 / pow(EL[0][i + 1], 2))

    # smoothing
    EL_N[N - 1] = EL[0][N - 1]
    VL_N[N - 1] = VL[0][N - 1]

    for i in range(0, N - 1) :
        i = N - 2 - i
        H = VL[0][i] / VL[1][i]

        EL_N[i] = EL[0][i] + H * (EL_N[i + 1] - EL[1][i])
        VL_N[i] = VL[0][i] + H * H * (VL_N[i + 1] - VL[1][i])
        COVL_N[i] = H * VL_N[i + 1]

    return [EL_N, VL_N, COVL_N]

####
# drawBRE
# draws the rate of event occurrence.
# arguments:
# spike_times: a spike train
# kalman: estimated rate
####

def drawBRE(spike_times, kalman) :
    xaxis = []
    yaxis = kalman[0][:]
    for i in range(0, len(spike_times) - 1) :
        xaxis.append((spike_times[i] + spike_times[i + 1]) / 2)

    plt.stackplot(xaxis, yaxis)
    plt.show()
