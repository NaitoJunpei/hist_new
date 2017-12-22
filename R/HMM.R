hmm <- function(spike_times) {
    ################################
    # determine the times of initiation and termination
    # bin size =  5*(inter-spike interval)
    ##############################

    max_value <- max(spike_times)
    min_value <- min(spike_times)
    onset <- min_value - 0.001 * (max_value - min_value)
    offset <- max_value + 0.001 * (max_value - min_value)
    bin_width <- (offset - onset) / length(spike_times) * 5

    ###############################
    # get_hmm_ratefunc: compute the transition between the hidden states
    # rate_hmm = (time, rate) given in a form of ndarray
    # drawHMM: draw a figure
    ##############################

    rate_hmm = getHmmRatefunc(spike_times, bin_width, max_value, min_value)

    drawHMM(spike_times, rate_hmm)
    return (rate_hmm)
}

drawHMM <- function(spike_times, rate_hmm) {
    xaxis <- rate_hmm[1, 1]
    yaxis <- rate_hmm[1, 2]
    tempx_old <- tempx <- rate_hmm[1, 1]
    tempy_old <- tempy <- rate_hmm[1, 2]

    for (i in 1:length(rate_hmm)) {
        tempx <- rate_hmm[i, 1]
        tempy <- rate_hmm[i, 2]
        if (tempy != tempx_old) {
            mid <- (tempx + tempx_old) / 2
            xaxis <- c(xaxis, mid, mid)
            yaxis <- c(yaxis, tempy_old, tempy)
        }
        tempx_old <- tempx_old
        tempy_old <- tempy
    }

    xaxis <- c(xaxis, rate_hmm[nrow(rate_hmm), 1])
    yaxis <- c(yaxis, rate_hmm[nrow(rate_hmm), 2])

    plot(xaxis, yaxis, type="l")
}

getHmmRatefunc <- function(spike_times, bin_width, max_value, min_value) {
    EMloop_sum <- 5000

    mat_A <- matrix(c(0.999, 0.001, 0.001, 0.999), nrow=2, ncol=2, byrow=TRUE)
    vec_pi <- c(0.5, 0.5)
    mean_rate <- length(spike_times) / (max_value - min_value)

    vec_lambda <- c(mean_rate * 0.75 * bin_width, mean_rate * 1.25 * bin_width)

    vec_spkt <- spike_times - min_value

    vec_Xi <- getVecXi(vec_spkt, bin_width)

    #########################################################
    # Optimizing the model parameters using the Baum-Welch algorithm
    #
    # updates parameters by hmm_E_step and hmm_M_step
    ########################################################

    res <- hmmEStep(vec_Xi, mat_A, vec_lambda, vec_pi)
    mat_Gamma <- res[[1]]
    mat_Xi <- res[[2]]

    mat_A_old <- mat_A
    vec_pi_old <- vec_pi
    vec_lambda_old <- vec_lambda

    #####################
    # Evaluation in the while loop
    # stops when the change in a model parameter becomes small
    # set flag=1 if the sum of change in a parameter sumcheck becomes smaller than some threshold
    # or stops if the loops are repeated so many times
    #####################

    loop <- 0
    flag <- 0
    while((loop <= EMloop_sum) && (flag == 0)) {
        res <- hmmMStep(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi)
        vec_pi_new <- res[[1]]
        vec_lambda_new <- res[[2]]
        mat_A_new <- res[[3]]

        vec_pi <- vec_pi_new
        vec_lambda <- vec_lambda_new
        mat_A <- mat_A_new

        sum_check <- 0.0
        num_state <- length(vec_pi)

        sum_check <- sum_check + sum(abs(vec_pi_old - vec_pi))
        sum_check <- sum_check + sum(abs(vec_lambda_old - vec_lambda))
        sum_check <- sum_check + sum(abs(mat_A_old - mat_A))

        if (sum_check / (1.0 * num_state * (num_state + 2)) < 1.0e-7) {
            flag <- 1
        }

        mat_A_old <- mat_A
        vec_pi_old <- vec_pi
        vec_lambda_old <- vec_lambda

        res <- hmmEStep(vec_Xi, mat_A, vec_lambda, vec_pi)
        mat_Gamma <- res[[1]]
        mat_Xi <- res[[2]]

        loop <- loop + 1
    }

    #############################################
    #
    #  Estimate an optimal sequence of states using the Viterbi algorithm
    # state is represented as 0 or 1 here
    #
    #############################################

    vec_hidden <- hmmViterbi(vec_Xi, mat_A, vec_lambda, vec_pi)

    ############################################
    # vec_hidden: 0 or 1, representing the hidden state
    # two states are transformed into the rates
    ############################################

    rate_func <- matrix(0, nrow=length(vec_Xi), 2)

    c_time = 0.0
    for (n in 1:length(vec_Xi)) {
        state_id <- vec_hidden[n]
        rate_func[n, 1] <- round(c_time * 100) / 100.0
        rate_func[n, 2] <- round(vec_lambda[floor(state_id) + 1] * 100) / (bin_width * 100.0)

        c_time <- c_time + bin_width
    }

    return (rate_func)
}

getVecXi <- function(vec_spkt, bin_width) {
    spkt_dura <- vec_spkt[length(vec_spkt)]
    bin_num <- ceiling(spkt_dura / bin_width)
    vec_Xi <- numeric(bin_num)

    ##########
    # counting spikes
    ##########
    for (x in vec_spkt) {
        bin_id <- floor(x / bin_width) + 1
        if (bin_id <= bin_num) {
            vec_Xi[bin_id] <- vec_Xi[bin_id] + 1
        }
    }

    return (vec_Xi)
}

hmmEStep <- function(vec_Xi, mat_A, vec_lambda, vec_pi) {
    mat_emission <- getMatEmission(vec_Xi, vec_lambda)
    res <- getAlphaC(mat_A, vec_pi, mat_emission)
    vec_C <- res[[1]]
    mat_alpha <- res[[2]]
    mat_beta <- getBeta(mat_A, vec_pi, mat_emission, vec_C)
    res <- getGammaXi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)
    mat_Gamma <- res[[1]]
    mat_Xi <- res[[2]]

    return (list(mat_Gamma, mat_Xi))
}

getMatEmission <- function(vec_Xi, vec_lambda) {
    mat_emission <- matrix(0, nrow=length(vec_Xi), ncol=length(vec_lambda))
    for (index_y in 1:length(vec_lambda)) {
        for (index_x in 1:length(vec_Xi)) {
            mat_emission[index_x, index_y] <- (vec_lambda[index_y])^(vec_Xi[index_x]) * exp(-1.0 * vec_lambda[index_y]) / factorial(vec_Xi[index_x])
        }
    }

    return (mat_emission)
}

getAlphaC <- function(mat_A, vec_pi, mat_emission) {
    num_of_states <- length(vec_pi)
    num_of_obs <- nrow(mat_emission)

    mat_alpha <- matrix(0, nrow=num_of_obs, ncol=num_of_states)
    coefficient <- numeric(num_of_obs)

    #####################
    # foward algorithm
    #
    # initialization
    mat_alpha[1,] <- vec_pi * mat_emission[1,]
    # innduction
    for (t in 1:(num_of_obs - 1)) {
        coefficient[t] <- sum(mat_alpha[t,])
        mat_alpha[t,] <- mat_alpha[t,] / coefficient[t]
        mat_alpha[t+1,] <- (mat_alpha[t,] %*% mat_A) * mat_emission[t+1,]
    }
    coefficient[num_of_obs] <- sum(mat_alpha[num_of_obs,])
    mat_alpha[num_of_obs,] <- mat_alpha[num_of_obs,] / coefficient[num_of_obs]
    return (list(coefficient, mat_alpha))
}

getBeta <- function(mat_A, vec_pi, mat_emission, vec_C) {
    num_of_states <- length(vec_pi)
    num_of_obs <- nrow(mat_emission)

    mat_beta <- matrix(0, nrow=num_of_obs, ncol=num_of_states)

    for (i in 1:num_of_states) {
        mat_beta[num_of_obs, i] <- 1.0
    }
    for (m in 2:num_of_obs) {
        
        n <- num_of_obs + 1 - m
        mat_beta_n1 <- mat_beta[n+1]
        mat_emission_n1 <- mat_emission[n + 1]
        mat_beta_n <- mat_beta[n]
        vec_C_n1 <- vec_C[n+1]
        for (i in 1:num_of_states) {
            mat_A_i <- mat_A[i,]
            sum_j <- sum(mat_beta_n1 * mat_emission_n1 * mat_A_i)
            mat_beta[n, i] <- (sum_j / vec_C_n1)
        }
    }

    return (mat_beta)
}

getGammaXi <- function(mat_A, mat_emission, mat_alpha, mat_beta, vec_C) {
    num_of_states <- ncol(mat_emission)
    num_of_obs <- nrow(mat_emission)
    mat_Gamma <- mat_alpha * mat_beta
    mat_Xi <- array(0, dim=c(num_of_obs - 1, num_of_states, num_of_states))
    for (t in 1:(num_of_obs - 1)) {
        # mat_Xi[t,] <- matrix(mat_alpha[t,], nrow=num_of_states, ncol=1) * mat_A * mat_emission[t + 1,] * mat_beta[t + 1,] / vec_C[t + 1]
        mat_Xi[t,,] <- apply(mat_A, MARGIN=2, FUN=function(x) { x * mat_alpha[t,]}) * mat_emission[t + 1,] * mat_beta[t + 1,] / vec_C[t + 1]
    }
    return (list(mat_Gamma, mat_Xi))
}

hmmMStep <- function(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi) {
    num_of_states <- nrow(mat_A)
    num_of_obs <- length(vec_Xi)

    vec_pi_new <- mat_Gamma[1,] / sum(mat_Gamma[1,])

    # vec_lambda_new <- rowSums(mat_Gamma * matrix(vec_Xi, nrow=num_of_obs, ncol=1)) / rowSums(mat_Gamma)
    vec_lambda_new <- colSums(apply(mat_Gamma, MARGIN=2, FUN=function(x) { x * vec_Xi })) / colSums(mat_Gamma)

    total_gamma_a = matrix(0, nrow=1, ncol=num_of_states)
    for (t in 1:(num_of_obs - 1)) {
        total_gamma_a <- total_gamma_a + mat_Gamma[t,]
    }

    total_xi <- colSums(mat_Xi)
    # mat_A_new <- total_xi / t(total_gamma_a)
    mat_A_new <- apply(total_xi, MARGIN=2, FUN=function(x) { x / total_gamma_a })

    return (list(vec_pi_new, vec_lambda_new, mat_A_new))
}

hmmViterbi <- function(vec_Xi, mat_A, vec_lambda, vec_pi) {
    mat_emission <- getMatEmission(vec_Xi, vec_lambda)
    num_of_states <- nrow(mat_A)
    num_of_obs <- length(vec_Xi)
    mat_hs_seq <- matrix(0, nrow=num_of_states, ncol=num_of_obs)
    vec_logp_seq <- numeric(num_of_states)

    for (j in 1:num_of_states) {
        mat_hs_seq[j, 1] <- j
        if (vec_pi[j] * mat_emission[1, j] == 0) {
            vec_logp_seq[j] <- -Inf
        } else {
            vec_logp_seq[j] <- log(vec_pi[j] * mat_emission[1, j], base=10)
        }
    }
    for (n in 1:num_of_obs) {
        mat_hs_seq_buf <- mat_hs_seq
        vec_logp_seq_buf <- vec_logp_seq

        for (j in 1:num_of_states) {
            vec_h_logprob_i <- numeric(num_of_states)
            for (i in 1:num_of_states) {
                vec_h_logprob_i[i] <- vec_logp_seq[i] + log(mat_emission[n, j] * mat_A[i][j], base=10)
            }

            max_element <- max(vec_h_logprob_i)
            max_pos <- which(vec_h_logprob_i == max_element)

            vec_logp_seq_buf[j] <- max_element
            mat_hs_seq_buf[j,] <- mat_hs_seq[max_pos,]
            print(mat_hs_seq_buf[j,])
            print(mat_hs_seq[max_pos,])
            print(max_pos)

            mat_hs_seq_buf[j, n] <- j
        }

        mat_hs_seq <- mat_hs_seq_buf
        vec_logp_seq <- vec_logp_seq_buf
    }

    max_element <- max(vec_logp_seq)
    max_pos <- which(vec_logp_seq == max_element)

    vec_hs_seq <- mat_hs_seq[max_pos,]

    return (vec_hs_seq)
}
