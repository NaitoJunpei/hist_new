##########
# HMM.py  returns the firing rate selected as an alternative hidden state.
#  needs libraries: (matplotlib, numpy, pandas). 

# Instruction
# put HMM.py in a folder on a path.
# import HMM
# then you may obtain HMM.().
# you need only HMM function.
# the function HMM take a spike train as an argument.
# spike train could be given by list or numpy.array.
# parameters are determined by the HMM and  a figure is drawn.
# references:
# Mochizuki and Shinomoto, Analog and digital codes in the brain
# Physical Review E (2014) 89:022705
# https://arxiv.org/abs/1311.4035
# Contact:
# Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp
##########
###########################################
# HMM_v2.py
# 2017/11/13 revised by Daisuke Endo
# endo.daisuke.63m@st.kyoto-u.ac.jp
###########################################

#############
# Instruction
#
# put HMM_v1.py in a folder on a path.
# import HMM_v1 as HMM
#  then executable as HMM.(function name)
#
# the function HMM takes a spike train as an argument
# spike train could be given by list or numpy.array
# read a text file of a spike train by data = np.loadtxt("data.txt")
# this program computes hidden variables given as instantaneous rate, using the HMM and draw a figure of the rate.
#
# you need libraries of matplotlib, numpy, and pandas
# references:
# Mochizuki and Shinomoto, Analog and digital codes in the brain
# https://arxiv.org/abs/1311.4035
# Contact:
# Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp
############

import matplotlib.pyplot as plt
import numpy as np
import math


def hmm(spike_times):
    """
    argument: spike_times: given in list or ndarray
    returns hidden states (rates) and draw the figure
    observed values is not binary (0,1) but the number of spikes in a given time bin.
    Poisson distribution is assumed for the number of spikes for a given time bin.
 
python HMM_v1.py:
    if __name__ == "__main__":
        data = np.loadtxt("data.txt")
        hmm(spike_times=data)

    # HMM.py: 
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import HMM_v1 as HMM
    data = np.loadtxt("data.txt")
    HMM.hmm(spike_times=data)
    """
    ################################
    # determine the times of initiation and termination
    # bin size =  5*(inter-spike interval)
    ##############################

    spike_times = np.array(list(spike_times))
    max_value = max(spike_times)
    min_value = min(spike_times)
    onset = min_value - 0.001 * (max_value - min_value)
    offset = max_value + 0.001 * (max_value - min_value)
    bin_width = (offset - onset) / len(spike_times) * 5

    ###############################
    # get_hmm_ratefunc: compute the transition between the hidden states
    # rate_hmm = (time, rate) given in a form of ndarray
    # drawHMM: draw a figure
    ##############################

    rate_hmm = get_hmm_ratefunc(spike_times, bin_width, max_value, min_value)

    drawHMM(spike_times, rate_hmm)

    return rate_hmm


def drawHMM(spike_times, rate_hmm):
    """
    draw a figure of the estimated state (given in a form of firing rate)
    arguments:
        spike_times : given in list or ndarray
        rate_hmm = (time, rate) determined by the Baum-Welch algorithm and the Viterbi algorithm in a form of ndarray
    returns:
        nothing, but draws a figure
    """
    ###################
    # output is shaped so that state changes appear vertical.
    ###################

    xaxis = [rate_hmm[0, 0]]
    yaxis = [rate_hmm[0, 1]]
    tempx_old = tempx = rate_hmm[0, 0]
    tempy_old = tempy = rate_hmm[0, 1]
    for i in range(0, len(rate_hmm) - 1):
        tempx, tempy = rate_hmm[i]
        if (tempy != tempy_old):
            mid = (tempx + tempx_old) / 2
            xaxis.append(mid)
            xaxis.append(mid)
            yaxis.append(tempy_old)
            yaxis.append(tempy)
        tempx_old = tempx
        tempy_old = tempy

    xaxis.append(rate_hmm[-1, 0])
    yaxis.append(rate_hmm[-1, 1])
    plt.stackplot(xaxis, yaxis)
    plt.xlim(xmin=min(xaxis), xmax=max(xaxis))
    plt.ylim(ymin=0)
    plt.show()


def get_hmm_ratefunc(spike_times, bin_width, max_value, min_value):
    """
    infers optimal states using the Baum-Welch algorithm and the Viterbi algorithm
   arguments:
        spike_times
        bin_width: bin size
        max_value: final spike time
        min_value: initial spike time
    returns:
        hidden states (rate_func). Here rate_func is given by a matrix of (time, state)
    example:
    rate_hmm = get_hmm_ratefunc(spike_times, bin_width, max_values, min_values)
    """
    #######################
    # set the initial values of the model parameters
    # vec_spkt: sets the initial spike time
    # vec_Xi: acquires the observed values
    # vec_Xi consists of the number of spikes (0, 1, 2, 3, ..) in each step.
    #######################
    EMloop_num = 5000

    mat_A = np.array([[0.999, 0.001], [0.001, 0.999]])
    vec_pi = np.array([0.5, 0.5])
    mean_rate = len(spike_times) / (max_value - min_value)

    vec_lambda = np.empty(2)

    vec_lambda[0] = (mean_rate * 0.75) * bin_width
    vec_lambda[1] = (mean_rate * 1.25) * bin_width

    vec_spkt = np.array([spike - min_value for spike in spike_times])

    vec_Xi = get_vec_Xi(vec_spkt, bin_width)

    #########################################################
    #
    # Optimizing the model parameters using the Baum-Welch algorithm
    #
    # updates parameters by hmm_E_step and hmm_M_step
    ########################################################

    mat_Gamma, mat_Xi = hmm_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)

    mat_A_old = mat_A
    vec_pi_old = vec_pi
    vec_lambda_old = vec_lambda

    #####################
    # Evaluation in the while loop
    # stops when the change in a model parameter becomes small
    # set flag=1 if the sum of change in a parameter sumcheck becomes smaller than some threshold
    # or stops if the loops are repeated so many times 
    #####################

    loop = 0
    flag = 0
    while(loop <= EMloop_num and flag == 0):
        vec_pi_new, vec_lambda_new, mat_A_new = hmm_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi)

        vec_pi = vec_pi_new
        vec_lambda = vec_lambda_new
        mat_A = mat_A_new

        sum_check = 0.0
        num_state = len(vec_pi)

        sum_check += sum(abs(vec_pi_old - vec_pi))

        sum_check += sum(abs(vec_lambda_old - vec_lambda))

        sum_check += sum(sum(abs(mat_A_old - mat_A)))

        if (sum_check / (1.0 * num_state * (num_state + 2)) < 1.0e-7):
            flag = 1
        mat_A_old = mat_A
        vec_pi_old = vec_pi
        vec_lambda_old = vec_lambda

        mat_Gamma, mat_Xi = hmm_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)

        loop += 1

    #############################################
    #
    #  Estimate an optimal sequence of states using the Viterbi algorithm
    # state is represented as 0 or 1 here 
    #
    #############################################

    vec_hidden = hmm_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi)

    ############################################
    # vec_hidden: 0 or 1, representing the hidden state
    # two states are transformed into the rates
    ############################################

    rate_func = np.empty([len(vec_Xi), 2])

    c_time = 0.0
    for n in range(0, len(vec_Xi)):
        state_id = vec_hidden[n]
        rate_func[n][0] = round(c_time * 100) / 100.0
        rate_func[n][1] = round(vec_lambda[int(state_id)] * 100) / (bin_width * 100.0)

        c_time += bin_width

    return rate_func


################################
#
# Function acquiring the observation sequence from a spike train
#
################################
def get_vec_Xi(vec_spkt, bin_width):
    """
    arguments:
        vec_spkt: spike time measured from the initial spike
        bin_width: bin size
    returns:
        vec_Xi: observation values consisting of spike counts in each bin.
    """
    spkt_dura = vec_spkt[len(vec_spkt) - 1]
    bin_num = int(math.ceil(spkt_dura / bin_width))
    vec_Xi = np.zeros(bin_num)

    ##########
    # counting spikes
    ##########
    for x in vec_spkt:
        bin_id = int(math.floor(x / bin_width))
        if (bin_id < bin_num):
            vec_Xi[bin_id] += 1
    return vec_Xi


###########################
#
# hmm_E_step
# computes expectation 
#
###########################
def hmm_E_step(vec_Xi, mat_A, vec_lambda, vec_pi):
    """
    arguments:
        vec_Xi: observation consisting of spike counts
        mat_A: transition matrix
        vec_lambda: spikes in each bin
        vec_pi: initial probabilities
    returns:
            mat_Gamma: a matrix consisting of P(state i at time t, vec_Xi |model) 
            mat_Xi: a matrix consisting of P(state i at time t and state j at time t+1, vec_Xi|model)
    """
    mat_emission = get_mat_emission(vec_Xi, vec_lambda)
    vec_C, mat_alpha = get_alpha_C(mat_A, vec_pi, mat_emission)
    mat_beta = get_beta(mat_A, vec_pi, mat_emission, vec_C)
    mat_Gamma, mat_Xi = get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)

    return mat_Gamma, mat_Xi


def get_mat_emission(vec_Xi, vec_lambda):
    """
    arguments:
        vec_Xi: observation
        vec_lambda: average spikes in each bin
    returns:
        mat_emission: probability of obtaining each observation in each state
        mat_emission: matrix consisting of (step,state)
    """
    ############
    # assuming the Poisson distribution
    # vec_lambda: parameter of the Poisson distribution representing the mean
    # gives the probability of having the observation
    #############
    mat_emission = np.array([[pow(x, y) * pow(math.e, -1.0 * x) / math.factorial(y) for x in vec_lambda] for y in vec_Xi])

    return mat_emission


def get_alpha_C(mat_A, vec_pi, mat_emission):
    """
    arguments:
        mat_A: transition matrix
        vec_pi: initial probability
        mat_emission: matrix consisting of the probability of having the observation (step,state)
    returns:
            mat_alpha: forward parameter
            vec_C(coefficient): scaling coefficient
    """

    num_of_states = len(vec_pi)
    num_of_obs = len(mat_emission)

    mat_alpha = np.empty((num_of_obs, num_of_states))
    coefficient = np.empty(num_of_obs)

    #####################
    # foward algorithm
    #
    # initialization
    mat_alpha[0, :] = vec_pi * mat_emission[0, :]
    # induction
    for t in range(0, num_of_obs - 1):  # do scaling
        coefficient[t] = np.sum(mat_alpha[t, :])
        mat_alpha[t, :] = mat_alpha[t, :] / coefficient[t]
        mat_alpha[t+1, :] = (mat_alpha[t, :] @ mat_A) * mat_emission[t+1, :]
    coefficient[num_of_obs - 1] = np.sum(mat_alpha[num_of_obs - 1, :])
    mat_alpha[num_of_obs - 1, :] = mat_alpha[num_of_obs - 1, :] / coefficient[num_of_obs - 1]

    return coefficient, mat_alpha


def get_beta(mat_A, vec_pi, mat_emission, vec_C):
    """
    arguments:
        mat_A: transition matrix
        vec_pi: initial probability
        mat_emission: matrix consisting of the probability of having the observation (step,state)
        vec_C: scaling coefficient obtained when computing alpha
    returns:
        mat_beta: backward parameter
    """
    ##############################
    # note given by Endo
    # here we do not introduce matrix calculation for the state
    # because the for-loop calculation was more rapid.
    # we leave here the code for matrix computation.
    ##############################
    ##############################
    # note:
    # when computing mat_beta[t, :], it is not divided by vec_C[t]
    # to avoid duplication when computing mat_Gamma
    # when computing mat_Xi, it was divided by vec_C[t+1] to make end meets. ##############################
    num_of_states = len(vec_pi)
    num_of_obs = len(mat_emission)

    mat_beta = np.zeros([num_of_obs, num_of_states])
    ###############
    # backward algorithm
    #
    # initialization
    for i in range(0, num_of_states):
        mat_beta[num_of_obs - 1][i] = 1.0
    # induction
    for m in range(1, num_of_obs):
        n = num_of_obs - 1 - m
        mat_beta_n1 = mat_beta[n + 1]
        mat_emission_n1 = mat_emission[n + 1]
        mat_beta_n = mat_beta[n]
        vec_C_n1 = vec_C[n + 1]
        for i in range(0, num_of_states):
            sum_j = 0.0
            mat_A_i = mat_A[i]
            for j in range(0, num_of_states):
                sum_j += mat_beta_n1[j] * mat_emission_n1[j] * mat_A_i[j]

            mat_beta_n[i] = (sum_j / vec_C_n1)

    return mat_beta
    '''
    # the following is the matrix computation
    num_of_states = len(vec_pi)
    num_of_obs = len(mat_emission)
    beta = np.empty((num_of_obs, num_of_states))
    ######################
    # backward algorithm
    # initialization beta
    beta[num_of_obs - 1, :] = 1
    # induction
    for t in range(num_of_obs - 1, 0, -1):
        beta[t - 1, :] = np.sum((mat_A * mat_emission[t, :] * beta[t, :]), axis=1) / vec_C[t]
    return beta'''


def get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C):
    """
      arguments:
        mat_A: transition matrix
        vec_pi: initial probability
        mat_emission: matrix consisting of the probability of having the observation (step,state)
        vec_C: scaling coefficient obtained when computing alpha
    returns:
            mat_Gamma: a matrix consisting of P(state i at time t, vec_Xi |model) 
            mat_Xi: a matrix consisting of P(state i at time t and state j at time t+1, vec_Xi|model)
    """
    num_of_states = len(mat_emission[0])
    num_of_obs = len(mat_emission)

    mat_Gamma = mat_alpha * mat_beta

    mat_Xi = np.empty((num_of_obs - 1, num_of_states, num_of_states))
    for t in range(0, num_of_obs - 1):
        mat_Xi[t, :, :] = mat_alpha[t, :].reshape(
            num_of_states, 1) * mat_A * mat_emission[t + 1, :] * mat_beta[t + 1, :] / vec_C[t+1]
    ################
    # note
    # procedure of dividing by vec_C when computing mat_Xi is complicated but
    # this is because mat_beta[t+1, :] was not divided by vec_C[t+1].
    ################
    return mat_Gamma, mat_Xi


###########################
#
# hmm_M_step
# this is a M step, maximizing the likelihood
#
###########################
def hmm_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi):
    """
    arguments:
        mat_A: transition matrix
        vec_pi: initial probability
        mat_Gamma: a matrix consisting of P(state i at time t, vec_Xi |model) 
        mat_Xi: a matrix consisting of P(state i at time t and state j at time t+1, vec_Xi|model)
    returns:
        updated parameters: vec_pi_new, vec_lambda_new, mat_A_new
    """
    num_of_states = len(mat_A)
    num_of_obs = len(vec_Xi)

    #############################
    #
    # computes new parameters
    #
    #############################
    #################
    # vec_pi
    ################
    vec_pi_new = mat_Gamma[0] / np.sum(mat_Gamma[0])

    #################
    # vec_lambda
    #################
    vec_lambda_new = np.sum(mat_Gamma * vec_Xi.reshape(num_of_obs, 1), axis=0) / np.sum(mat_Gamma, axis=0)

    ###############
    # mat_A
    ###############
    total_gamma_a = np.zeros((1, num_of_states))
    for t in range(0, num_of_obs - 1):
        total_gamma_a += mat_Gamma[t, :]
    total_xi = np.sum(mat_Xi, axis=0)  # (i,j): 2 by 2 matrix
    mat_A_new = total_xi / total_gamma_a.T  # division between i-th ones 
    return vec_pi_new, vec_lambda_new, mat_A_new


#####################################
#
#  a function for determining an optimal state sequence with the Viterbi algorithm
#
#####################################
def hmm_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi):
    """
    arguments:
        vex_Xi: observation sequence
        mat_A: transition matrix
        vec_lambda: average spike rate in each bin
        vec_pi: initial probability
    returns:
        vec_hs_seq: optimal state sequence 
    """
    mat_emission = get_mat_emission(vec_Xi, vec_lambda)
    num_of_states = len(mat_A)
    num_of_obs = len(vec_Xi)
    mat_hs_seq = np.zeros([num_of_states, num_of_obs])
    vec_logp_seq = np.zeros(num_of_states)

    for j in range(0, num_of_states):
        mat_hs_seq[j][0] = j
        if (vec_pi[j] * mat_emission[0][j] == 0):
            vec_logp_seq[j] = -np.inf
        else:
            vec_logp_seq[j] = math.log(vec_pi[j] * mat_emission[0][j]) / math.log(10)

    for n in range(1, num_of_obs):
        # copy the seq. up to n - 1
        mat_hs_seq_buf = mat_hs_seq.copy()
        vec_logp_seq_buf = vec_logp_seq.copy()

        for j in range(0, num_of_states):
            vec_h_logprob_i = np.zeros(num_of_states)
            for i in range(0, num_of_states):
                vec_h_logprob_i[i] = vec_logp_seq[i] + math.log(mat_emission[n][j] * mat_A[i][j]) / math.log(10)

            max_element = max(vec_h_logprob_i)
            max_pos = np.where(vec_h_logprob_i == max_element)[0][0]

            vec_logp_seq_buf[j] = max_element
            mat_hs_seq_buf[j] = mat_hs_seq[max_pos].copy()

            mat_hs_seq_buf[j][n] = j

        mat_hs_seq = mat_hs_seq_buf.copy()
        vec_logp_seq = vec_logp_seq_buf.copy()

    max_element = max(vec_logp_seq)

    max_pos = np.where(vec_logp_seq == max_element)[0][0]

    vec_hs_seq = mat_hs_seq[max_pos].copy()

    return vec_hs_seq
