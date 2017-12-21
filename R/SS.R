SS <- function(spike_times) {
    max_value <- max(spike_times)
    min_value <- min(spike_times)
    onset <- min_value - 0.001 * (max_value - min_value)
    offset <- max_value + 0.001 * (max_value - min_value)

    # computes the cost function by changing the number of bins
    # adopts the number of bins that minimizes the cost function

    for (bin_num in 1:500) {
        cost <- costAv(spike_times, onset, offset, bin_num, 10)
        if (bin_num == 1 || cost < cost_min) {
            cost_min = cost
            optimal_bin_num = bin_num
        }
    }

    drawSS(spike_times, optimal_bin_num, onset, offset)
    return (optimal_bin_num)
}

########## 
# costF 
# computes the cost function defined by Shimazaki and Shinomoto

# arguments:
# spike_times: spike train
# start: time of the initial spike
# end: time of the final spike
# bin_num: number of bins

# returns the cost function 
########## 

costF <- function(spike_times, start, end, bin_num) {
    bin_width <- (end - start) / bin_num
    histogram <- hist(x=spike_times, plot=FALSE, breaks=seq(from=start, to=end, bin_width))
    histogram <- histogram$counts
    av <- mean(histogram)
    va <- mean(histogram * histogram)

    return ((2.0 * av - (va - av * av)) / (bin_width * bin_width))
}

########## 
# costAv
# computes an average cost function with respect to initial binning positions.

# arguments:
# spike_times: spike train
# onset: time of an initial spike
# offset: time of a final spike
# bin_num: the number of bins
# times: the number of initial binning positions

# returns the averaged cost function
##########

costAv <- function(spike_times, onset, offset, bin_num, times) {
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
        temp <- temp + costF(TT_cat, start, end, bin_num)
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

drawSS <- function(spike_times, optimal_bin_num, start, end) {
    bin_width <- (end - start) / optimal_bin_num
    hist(x=spike_times, breaks=seq(from=start, to=end, bin_width))
}

