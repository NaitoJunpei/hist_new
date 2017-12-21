SS <- function(spike_time)
{
min_sp <- min(spike_time)   # the time of the first spike
max_sp <- max(spike_time)   # the time of the last spike
onset = spike_time[1] - 0.001 * (spike_time[length(spike_time)] - spike_time[1])
offset = spike_time[length(spike_time)] + 0.001 * (spike_time[length(spike_time)] - spike_time[1])

N <- 1:50                  # num of bins vector
D <- (max_sp-min_sp)/N      # bin size  vector

Cost <- numeric(length(N))     # cost function vector 

SN <- 10                    # num of partitioning positions for shift average

Lv <- Calc_Lv( spike_time ) # Lv of the spike train


# Computing the cost function

for (i in N){
  shift <- seq(-D[i]/2,D[i]/2,length=SN)
  w <- numeric(SN)
  for (p in 1 : SN){
    start = onset + shift[i]
    end = offset + shift[i]
    count <- numeric(i)
    for (spike in spike_time[spike_time>start & spike_time<end]) {
      count[floor((spike-shift[p]) / D[i])+1] <- count[floor((spike-shift[p]) / D[i])+1] +1
    }
    # num of spikes in each bin (Step 1 in ref[1] pp 3129-3130)


    # cost function for num of bins = N(i) (Step 3,4 in ref[1] pp 3129-3130)
    w[p] <- mean(2.0*count)-var(count)
  }
  Cost[i] <- sum(w)/D[i]^2
}
  # Minimizing the cost function (Step 5 in ref[1] pp 3129-3130)
  OptN <- which.min(Cost)
  hist(spike_time,breaks=OptN)

return(D[OptN])
}
