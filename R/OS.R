##########
# OS.R computes the optimal number of bins of time-histogram based on the optimization method proposed by Omi and Shinomoto, which may be applicable to non-Poisson spike trains. 

# Instruction
# you need only OS function.
# the function OS takes a spike train as an argument. 
# spike train could be given by a list or numpy.array.
# the program selects the optimal bin size for a given spike train and draws the histogram. 
# references:
#  Takahiro Omi & Shigeru Shinomoto, "Optimizing time histograms for non-Poissonian spike trains", Neural Computation 23, 3125 (2011).
# Contact:
# Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp
##########
# OS.R
# Junpei Naito 2017/12/21
##########


OS <- function(spike_times) {
    max_value <- max(spike_times)
    min_value <- min(spike_times)
    onset <- min_value - 0.001 * (max_value - min_value)
    offset <- max_value + 0.001 * (max_value - min_value)
    lv <- 0
    ISI <- diff(spike_times)

    # computes the firing irregularity Lv
    for (i in 1:(length(ISI) - 1)) {
        interval1 <- ISI[i]
        interval2 <- ISI[i + 1]

        if (interval1 + interval2 != 0) {
            lv <- lv + 3 * (interval1 - interval2)^2 / ((interval1 + interval2)^2 * (length(spike_times) - 2))
        } else {
            lv <- lv + 3 / (length(spike_times) - 2)
        }
    }

    # computes the cost function by changing the number of bins
    # adopts the number of bins that minimizes the cost function

    for (bin_num in 1:500) {
        cost <- costAv(spike_times, onset, offset, lv, bin_num, 10)
        if (bin_num == 1 || cost < cost_min) {
            cost_min = cost
            optimal_bin_num = bin_num
        }
    }

    drawOS(spike_times, optimal_bin_num, onset, offset)
    return (optimal_bin_num)
}

########## 
# costF 
# computes the cost function defined by Shimazaki and Shinomoto

# arguments:
# spike_times: spike train
# start: time of the initial spike
# end: time of the final spike
# lv: the value of local variation Lv, which measures the spiking irregularity
# bin_num: number of bins

# returns the cost function 
########## 

costF <- function(spike_times, start, end, lv, bin_num) {
    bin_width <- (end - start) / bin_num
    histogram <- hist(x=spike_times, plot=FALSE, breaks=seq(from=start, to=end, bin_width))
    histogram <- histogram$counts

    fano <- 2.0 * lv / (3.0 - lv)

    av <- mean(histogram)
    va <- mean(histogram * histogram)
    w_av <- mean(histogram * fano)
    fano_bin <- histogram
    fano_bin[fano_bin <= 2] <- 1.0
    fano_bin[fano_bin > 2] <- fano
    
    return ((2.0 * mean(histogram * fano_bin) - (va - av * av)) / (bin_width * bin_width))
}

########## 
# costAv
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

costAv <- function(spike_times, onset, offset, lv, bin_num, times) {
    temp <- 0.0
    bin_width <- (offset - onset) / bin_num
    TT = c(spike_times, (spike_times + (offset - onset)))

    # averages the cost with respect to the starting positions.
    # times: number of starting positions.

    for (i in 0:(times - 1)) {
        start <- onset + i * bin_width / times
        end <- offset + i * bin_width / times
        TT_cat <- TT[start < TT]
        TT_cat <- TT_cat[TT_cat < end]
        temp <- temp + costF(TT_cat, start, end, lv, bin_num)
    }
    
    return (temp / times)
}

##########
# drawSS
# draws a histogram

# arguments:
# spike_times: a spike train
# optimal_bin_num: an optimal number of bins
########## 

drawOS <- function(spike_times, optimal_bin_num, start, end) {
    bin_width <- (end - start) / optimal_bin_num
    hist(x=spike_times, breaks=seq(from=start, to=end, bin_width))
}

