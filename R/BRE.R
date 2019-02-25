##########
# BRE.R returns Bayesian estimation of the rate of event occurrence.

# Instruction
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
# BRE.R
# Junpei Naito 2017/12/27
##########

BRE <- function(spike_times) {
    sort(spike_time)    # sort
    max_value <- max(spike_times)
    min_value <- min(spike_times)

    ISI <- diff(spike_times)
    mu <- length(spike_times) / (max_value - min_value)
    beta0 <- mu^(-3)
    beta <- EMMethod(ISI, beta0)
    kalman <- kalmanFilter(ISI, beta)

    drawBRE(spike_times, kalman)
}

####
# EMmethod
# estimates a parameter beta with the EM algorithm.

# arguments:
# ISI: inter spike interval
# beta0: initial value of the parameter beta

# returns the estimated beta.
####

EMMethod <- function(ISI, beta0) {
    N <- length(ISI)
    beta <- 0
    beta_new <- beta0

    for (j in 1:100) {
        beta <- beta_new
        kalman <- kalmanFilter(ISI, beta)

        beta_new <- 0
        t0 <- 0

        indexes <- which(ISI[1:(N - 1)] > 0)

        beta_new <- sum((kalman[2,indexes + 1] + kalman[2,indexes] - 2 * kalman[3,indexes]
            + (kalman[1,indexes + 1] - kalman[1,indexes])
            * (kalman[1,indexes + 1] - kalman[1,indexes])) / ISI[indexes])

        t0 <- N - 1 - length(indexes)
        beta_new <- (N - t0 - 1) / (2 * beta_new)
    }
    
    return (beta_new)

}

####
# KalmanFilter
# estimates the rate of event occurrence with the Kalman filtering.
# arguments:
# ISI: inter spike interval
# beta: the parameter
# returns the rate of event occurrence.
####

kalmanFilter <- function(ISI, beta) {
    N <- length(ISI)
    IEL <- N / sum(ISI)
    IVL <- (IEL / 3)^2
    A = IEL - ISI[1] * IVL
    EL <- matrix(0, nrow=2, ncol=N)
    VL <- matrix(0, nrow=2, ncol=N)

    EL_N <- numeric(N)
    VL_N <- numeric(N)
    COVL_N <- numeric(N)

    EL[1,1] <- (A + sqrt(A * A + 4 * IVL)) / 2
    VL[1,1] <- 1 / (1 / IVL + 1 / (EL[1,1])^2)

    # prediction and filtering
    for (i in 1:(N - 1)) {
        EL[2,i] <- EL[1,i]
        VL[2,i] <- VL[1,i] + ISI[i] / (2 * beta)

        A <- EL[2,i] - ISI[i + 1] * VL[2,i]
        EL[1,i + 1] <- (A + sqrt(A * A + 4 * VL[2,i])) / 2
        VL[1,i + 1] <- 1 / (1 / VL[2,i] + 1 / (EL[1,i + 1])^2)
    }

    # smoothing
    EL_N[N] <- EL[1,N]
    VL_N[N] <- VL[1,N]

    for (i in 1:(N - 1)) {
        index <- N - i
        H <- VL[1,index] / VL[2,index]

        EL_N[index] <- EL[1,index] + H * (EL_N[index + 1] - EL[2,index])
        VL_N[index] <- VL[1,index] + H * H * (VL_N[index + 1] - VL[2,index])
    }

    COVL_N = c((VL[1,1:(N - 1)] / VL[2,1:(N-1)]) * VL_N[2:N], 0)

    return (matrix(c(EL_N, VL_N, COVL_N), nrow=3, ncol=N, byrow=TRUE))
}

drawBRE <- function(spike_times, kalman) {
    xaxis <- c()
    yaxis <- kalman[1,]
    
    for (i in 1:(length(spike_times)- 1)) {
        xaxis <- c(xaxis, (spike_times[i] + spike_times[i + 1]) / 2)
    }

    plot(xaxis, yaxis, type="l", ylim=range(0,yaxis))
}
