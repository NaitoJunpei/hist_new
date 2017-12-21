HMM <- function(spike_time)
{
# Function `HMM_v3' returns the firing rate selected as an alternative hidden state.
# Original paper:
# Mochizuki and Shinomoto, Analog and digital codes in the brain
# https://arxiv.org/abs/1311.4035
#
# Input argument
# spike_time:    Sample data vector.
#
# Output argument
# rate_func: ratefunction.
#            2D array stores
#            1: begining time of each bins in second
#            2: rate of each bin
# made by Kazuki Nakamura 20171221
# Contact: Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp

# determine the times of initiation and termination
# bin size =  5*(inter-spike interval)
onset <- spike_time[1] - 0.001 * (spike_time[length(spike_time)] - spike_time[1])
offset <- spike_time[length(spike_time)] + 0.001 * (spike_time[length(spike_time)] - spike_time[1])
optw <- (offset-onset)/(length(spike_time)) * 5

# input: sample data vector and bin size
# compute the rate in each bin
rate_func <- get_hmm_ratefunc(spike_time, optw)
return(rate_func)

# draw a graph
drawHMM(rate_func)
}

# Function acquiring the observation sequence from a spike train
#
# arguments:
#   vec_spkt: spike time measured from the initial spike
#   bin_width: bin size
# returns:
#   vec_Xi: observation values consisting of spike counts in each bin.
get_vec_Xi <- function(vec_spkt, bin_width){
  bin_num <- ceiling(vec_spkt[length(vec_spkt)]/bin_width)
  vec_Xi <- numeric(bin_num)
  # counting spikes
  for(i in 1:length(vec_spkt)){
    bin_id <- trunc(vec_spkt[i]/bin_width)+1
    if(bin_id<bin_num+1){
      vec_Xi[bin_id] <- vec_Xi[bin_id]+1
    }
  }
  return(vec_Xi)
}

# Function to get emisson probs.
# assuming the Poisson distribution
# vec_lambda: parameter of the Poisson distribution representing the mean
# gives the probability of having the observation
#
# arguments:
#   vec_Xi: observation
#   vec_lambda: average spikes in each bin
# returns:
#   mat_emission: probability of obtaining each observation in each state
#   mat_emission: matrix consisting of (step,state)
get_mat_emission <- function(vec_Xi, vec_lambda){
  mat_emission <- matrix(numeric(length(vec_Xi)*length(vec_lambda)), nrow=length(vec_Xi), ncol=length(vec_lambda))
  for(n in 1:length(vec_Xi)){
    for(i in 1:length(vec_lambda)){
      mat_emission[n,i] <- vec_lambda[i]^vec_Xi[n]*exp(-1.0*vec_lambda[i]/factorial(vec_Xi[n]))
    }
  }
  return(mat_emission)
}

# Function to get alpha and C
#
# arguments:
#   mat_A: transition matrix
#   vec_pi: initial probability
#   mat_emission: matrix consisting of the probability of having the observation (step,state)
# returns:
#   mat_alpha: forward parameter
#   vec_C(coefficient): scaling coefficient
get_alpha_C <- function(mat_A, vec_pi, mat_emission)
{
  num_of_states <- length(vec_pi)
  num_of_obs <- dim(mat_emission)[1]
  
  vec_C <- numeric(num_of_obs)
  mat_alpha <- matrix(numeric(num_of_obs*num_of_states),nrow=num_of_obs, ncol = num_of_states)
  # n=1
  mat_alpha[1,] <- t(mat_emission[1,]) * vec_pi
  vec_C[1] <- sum(mat_alpha[1,])
  mat_alpha[1,] <- mat_alpha[1,]/vec_C[1]
  # n>1
  for(n in 2:num_of_obs){
    for(i in 1:num_of_states){
      sum_j <- sum(t(mat_alpha[n-1,]) * mat_A[,i])
      mat_alpha[n,i] <- mat_emission[n,i]*sum_j
    }
    vec_C[n] <- sum(mat_alpha[n,])
    mat_alpha[n,] <- mat_alpha[n,]/vec_C[n]
  }
  return(list(vec_C, mat_alpha))
}

# Function to get beta
#
# arguments:
#   mat_A: transition matrix
#   vec_pi: initial probability
#   mat_emission: matrix consisting of the probability of having the observation (step,state)
#   vec_C: scaling coefficient obtained when computing alpha
# returns:
#   mat_beta: backward parameter
get_beta <- function(mat_A, vec_pi, mat_emission, vec_C)
{
  num_of_states <- length(vec_pi)
  num_of_obs <- dim(mat_emission)[1]
  
  # initialize
  mat_beta <- matrix(numeric(num_of_obs*num_of_states),nrow=num_of_obs,ncol=num_of_states)
  # n=N
  mat_beta[num_of_obs,]=1.0
  # n<N
  for(m in 2:num_of_obs){
    n <- num_of_obs+1-m
    for(i in 1:num_of_states){
      sum_j <- sum(mat_beta[n+1,] * mat_emission[n+1,] * mat_A[i,])
      mat_beta[n,i] <- sum_j/vec_C[n+1]
    }
  }
  return(mat_beta)
}

# Function to get Gamma
#
# arguments:
#   mat_A: transition matrix
#   vec_pi: initial probability
#   mat_emission: matrix consisting of the probability of having the observation (step,state)
#   vec_C: scaling coefficient obtained when computing alpha
# returns:
#   mat_Gamma: a matrix consisting of P(state i at time t, vec_Xi |model)
#   mat_Xi: a matrix consisting of P(state i at time t and state j at time t+1, vec_Xi|model)
get_Gamma_Xi <- function(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)
{
  num_of_obs <- dim(mat_emission)[1]
  num_of_states <- dim(mat_emission)[2]
  # gamma matrix
  mat_Gamma <- mat_alpha * mat_beta
  mat_Xi <- array(numeric((num_of_obs-1)*num_of_states*num_of_states),dim=c(num_of_obs-1, num_of_states, num_of_states))
  for(m in 1:(num_of_obs-1)){
    mat_Xi[m,,] <- mat_A %*% t(mat_alpha[m,,drop=F]) * mat_emission[m+1,] * mat_beta[m+1,] / vec_C[m+1]
  }
  return(list(mat_Gamma, mat_Xi))
}

# HMM_E_step computes expectation
#
# arguments:
#   vec_Xi: observation consisting of spike counts
#   mat_A: transition matrix
#   vec_lambda: spikes in each bin
#   vec_pi: initial probabilities
# returns:
#   mat_Gamma: a matrix consisting of P(state i at time t, vec_Xi |model)
#mat_Xi: a matrix consisting of P(state i at time t and state j at time t+1, vec_Xi|model)
HMM_E_step <- function(vec_Xi, mat_A, vec_lambda, vec_pi)
{
  mat_emission <- get_mat_emission(vec_Xi, vec_lambda)
  alpha_C <- get_alpha_C(mat_A, vec_pi, mat_emission)
  vec_C <- alpha_C[[1]]
  mat_alpha <- alpha_C[[2]]
  mat_beta <- get_beta(mat_A, vec_pi, mat_emission, vec_C)
  Gamma_Xi <- get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)
  mat_Gamma <- Gamma_Xi[[1]]
  mat_Xi <- Gamma_Xi[[2]]
  return(list(mat_Gamma,mat_Xi))
}

# HMM_M_step maximize the likelihood
#
# arguments:
#   mat_A: transition matrix
#   vec_pi: initial probability
#   mat_Gamma: a matrix consisting of P(state i at time t, vec_Xi |model)
#   mat_Xi: a matrix consisting of P(state i at time t and state j at time t+1, vec_Xi|model)
# returns:
#   updated parameters: vec_pi_new, vec_lambda_new, mat_A_new
HMM_M_step <- function(vec_Xi,  mat_A,  vec_lambda,  vec_pi,  mat_Gamma,  mat_Xi)
{
  num_of_states <- dim(mat_A)[1]
  num_of_obs <- length(vec_Xi)
  
  # maximize wrt pi vector
  pi_denom <- sum(mat_Gamma[1,])
  vec_pi_new <- mat_Gamma[1,]/pi_denom
  
  # maximize wrt lambda vector
  vec_lambda_new <- numeric(num_of_states)
  for(k in 1:num_of_states){
    lambda_denom <- sum(mat_Gamma[,k])
    lambda_nume <- sum(mat_Gamma[,k]*vec_Xi)
    
    if(lambda_denom==0.0){
      vec_laombda_new[k]=0.0
    }else{
      vec_lambda_new[k] <- lambda_nume/lambda_denom
    }
  }
  #maximize wrt A matrix
  mat_A_new <- matrix(numeric(num_of_states*num_of_states),ncol=num_of_states,nrow=num_of_states)
  for(j in 1:num_of_states){
    A_denom <- sum(sum(mat_Xi[,j,]))
    for(k in 1:num_of_states){
      A_nume <- 0.0
      for(n in 1:(num_of_obs-1)){
        A_nume <- A_nume+mat_Xi[n,j,k]
        if(A_denom==0){
          mat_A_new[j,k] <- 0.0
        }else{
          mat_A_new[j,k] <- A_nume/A_denom
        }
      }
    }
  }
  return(list(vec_pi_new, vec_lambda_new, mat_A_new))
}

# HMM_Viterbi
#   function for determining an optimal state squence with the Viterbi algorithm
#
# arguments:
#   vex_Xi: observation sequence
#   mat_A: transition matrix
#   vec_lambda: average spike rate in each bin
#   vec_pi: initial probability
# returns:
#   vec_hs_seq: optimal state sequence
HMM_Viterbi <- function(vec_Xi, mat_A, vec_lambda, vec_pi)
{
  mat_emission <- get_mat_emission(vec_Xi, vec_lambda)
  
  num_of_states <- dim(mat_A)[1]
  num_of_obs <- length(vec_Xi)
  
  mat_hs_seq <- matrix(numeric(num_of_states*num_of_obs),nrow=num_of_states, ncol=num_of_obs)
  vec_logp_seq <- numeric(num_of_states)
  
  # n=1
  for(j in 1:num_of_states){
    mat_hs_seq[j,1] <- j
    vec_logp_seq[j] <- log(vec_pi[j]*mat_emission[1,j])/log(10)
  }
  # n>1
  for(n in 2:num_of_obs){
    # copy the seq.
    mat_hs_seq_buf <- mat_hs_seq
    vec_logp_seq_buf <- vec_logp_seq
    # nth node -> j
    for(j in 1:num_of_states){
      # n-1th node -> j
      # compute logp for i->j trans
      vec_h_logprob_i <- numeric(num_of_states)
      vec_h_logprob_i <- vec_logp_seq+log(mat_emission[n,j]*mat_A[,j])/log(10)
      
      # get max logp
      max_element <- max(vec_h_logprob_i)
      max_pos <- which.max(vec_h_logprob_i)
      vec_logp_seq_buf[j] <- max_element
      mat_hs_seq_buf[j,] <- mat_hs_seq[max_pos,]
      mat_hs_seq_buf[j,n] <- j
    }
    # updata the seq
    mat_hs_seq <- mat_hs_seq_buf
    vec_logp_seq <- vec_logp_seq_buf
  }
  max_element <- max(vec_logp_seq)
  max_pos <- which.max(vec_logp_seq)
  vec_hs_seq <- mat_hs_seq[max_pos,]
  return(vec_hs_seq)
}

# get_hmm_ratefunc ompute the transition between the hidden states
#
# infers optimal states using the Baum-Welch algorithm and the Viterbi algorithm
# arguments:
#   spike_times
#   bin_width: bin size
# returns:
#   hidden states (rate_func). Here rate_func is given by a matrix of (time, state)
#   2D array stores
#   1: begining time of each bins in second
#   2: rate of each bin
get_hmm_ratefunc <- function(spike_time, bin_width)
{
# set the initial values of the model parameters
# vec_spkt: sets the initial spike time
# vec_Xi: acquires the observed values
# vec_Xi consists of the number of spikes (0, 1, 2, 3, ..) in each step.
  EMloop_num <- 5000    # number of EM iteration
  mat_A <- matrix(c(0.999, 0.001, 0.001, 0.999), ncol=2, nrow=2)
  vec_pi <- c(0.5, 0.5)
  mean_rate <- length(spike_time)/(spike_time[length(spike_time)]-spike_time[1])
  vec_lambda <- c((mean_rate*0.75)*bin_width, (mean_rate*1.25)*bin_width)
  
  vec_spkt <- numeric(length(spike_time))
  for(i in 1:length(spike_time)){
    vec_spkt[i] <- spike_time[i] - spike_time[1]
  }
  vec_Xi <- get_vec_Xi(vec_spkt, bin_width)
  
  # Optimizing the model parameters using the Baum-Welch algorithm
  # updates parameters by hmm_E_step and hmm_M_step
  Gamma_Xi <- HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)
  mat_Gamma <- Gamma_Xi[[1]]
  mat_Xi <- Gamma_Xi[[2]]
  mat_A_old <- mat_A
  vec_pi_old <- vec_pi
  vec_lambda_old <- vec_lambda
  # Evaluation in the while loop
  # stops when the change in a model parameter becomes small
  # set flag=1 if the sum of change in a parameter sumcheck becomes smaller than some threshold
  # or stops if the loops are repeated so many times
  loop <- 0
  flag <- 0
  while(loop<=EMloop_num & flag==0){
    result <- HMM_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi)
    vec_pi_new <- result[[1]]
    vec_lambda_new <- result[[2]]
    mat_A_new <- result[[3]]
    
    vec_pi <- vec_pi_new
    vec_lambda <- vec_lambda_new
    mat_A <- mat_A_new
    
    sum_check <- 0.0
    num_state <- length(vec_pi)
    sum_check <- sum(sum(abs(mat_A_old-mat_A)))
    sum_check <- sum_check+sum(abs(vec_pi_old-vec_pi))
    sum_check <- sum_check+sum(abs(vec_lambda_old-vec_lambda))
    
    if(sum_check/(1.0*num_state*(num_state+2))<10^(-7)){
      flag <- flag+1
    }
    mat_A_old <- mat_A
    vec_pi_old <- vec_pi
    vec_lambda_old <- vec_lambda
    
    Gamma_Xi <- HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)
    mat_Gamma <- Gamma_Xi[[1]]
    mat_Xi <- Gamma_Xi[[2]]
    
    loop <- loop+1
  }
  #  Estimate an optimal sequence of states using the Viterbi algorithm
  # state is represented as 0 or 1 here
  vec_hidden <- numeric(length(vec_Xi))
  vec_hidden <- HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi)
  
  # vec_hidden: 0 or 1, representing the hidden state
  # two states are transformed into the rates
  rate_func <- matrix(numeric(length(vec_Xi)*2),ncol=2,nrow=length(vec_Xi))
  onset <- spike_time[1] - 0.001 * (spike_time[length(spike_time)] - spike_time[1])
  c_time <- onset
  for(n in 1:length(vec_Xi)){
    state_id <- vec_hidden[n]
    rate_func[n,1] <- round(c_time*100)/100.0
    rate_func[n,2] <- round(vec_lambda[state_id]*100)/(bin_width*100.0)
    c_time <- c_time+bin_width
  }
  return(rate_func)
}


# draw a figure of the estimated state (given in a form of firing rate)
# arguments:
#   spike_times : given in list or ndarray
#   rate_hmm = (time, rate) determined by the Baum-Welch algorithm and the Viterbi algorithm in a form of ndarray
# returns:
#   nothing, but draws a figure
drawHMM <- function(rate_func){
  
}
              