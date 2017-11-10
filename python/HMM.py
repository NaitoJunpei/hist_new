###########################################
#
# 元のpythonの改訂版
#
# forループを行列計算できる所は変更してます
# 元のコードはコメントアウトしてあります
# 計算手順は買かえてないです。
# 改訂する時の都合でEMloop_num = 1000 に変更してあります。プログラムにミスがあると最後まで回って時間がかかるためです。
# 細かいことですが、=の位置を揃えるのはpythonのコーディングの規約に反してます
# 元のコメントは残してますが、追加で新しいコメントを入れてあります
#
# 11/9改訂
# 改訂者：遠藤大輔
###########################################

######################
# 以下から元のコードに少し手を加えたもの
######################


##########
# HMM.pyを実行するには、matplotlib、numpy、pandasライブラリが必要です

# 使い方
# HMM.pyを、パスが通っているフォルダに置き、
# import HMM
# をすると、ファイル内の関数が、HMM.(関数名)の形で実行可能になります。

# ユーザーが使用するのはHMM関数のみで十分です。
# HMM関数は、spike列を引数に取ります。
# spike列の形式は、list、numpy.arrayなどが利用可能です。
# 隠れマルコフモデルを使ってパラメータを推測し、グラフを描画します。
# 値を返しません。
##########

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import time

def HMM(spike_times):
    """
    spike_times=spike時刻をlist,ndarrayなどに入れたもの
    結果：
    隠れstateの遷移結果の値とそのグラフを出力
    計算にかかった時刻も出力
    """
    ################################
    # スタート時刻、終了時刻を決める
    # bin幅は1つのbinに平均して5個のspikeが入るように決めている
    ##############################
    t1 = time.time()
    spike_times = np.array(list(spike_times))
    max_value   = max(spike_times)
    min_value   = min(spike_times)
    onset       = min_value - 0.001 * (max_value - min_value)
    offset      = max_value + 0.001 * (max_value - min_value)
    bin_width   = (offset - onset) / len(spike_times) * 5

    ###############################
    #
    # get_hmm_ratefuncでHMMの遷移の様子を求めます。
    # rate_hmm = (time, rate)の形のndarray
    # drawHMMで実際に描画します。
    #
    ##############################

    rate_hmm = get_hmm_ratefunc(spike_times, bin_width, max_value, min_value)
    t2 = time.time()
    print("baumwelch+viterbi time",t2-t1)
    drawHMM(spike_times, rate_hmm)
    return rate_hmm

####
# 隠れマルコフモデルで推定した値の描画を行います。

# 引数
# spike_times: スパイク列
# rate_hmm: 隠れマルコフモデルで推定した値

# 返り値
# なし
####


def drawHMM(spike_times, rate_hmm):
    """
    HMMで推定した値の描画を行う
    rate_hmm = バウムウェルチとビタビで推定された値。(time, rate)の形のndarray
    結果：
    グラフを出力
    """
    ###################
    # stateの変化が斜め線にならないように、いくらか補正しています。
    # forのなかがその操作です。
    # ややこしいですが、見栄えを良くするためです
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
    spike_times=spike時刻の列
    bin_width=bin幅
    max_value=spike時刻の最後の値
    min_value=spike時刻の最初の値
    結果：
    隠れstateの変化(rate_func)を出力
    rate_func = (time, state)の行列
    例：
    rate_hmm = get_hmm_ratefunc(spike_times, bin_width, max_values, min_values)
    """
    #######################
    # モデルパラメータの初期値を設定しています。
    # vec_spktでspike時刻列を初めのspike時刻がスタートになるよう補正しています。
    # vec_Xiでは観測列を得ています。
    # vec_Xiには各stepでのspikeの個数が入っています。
    #######################
    EMloop_num = 1000

    mat_A      = np.array([[0.999, 0.001], [0.001, 0.999]])
    vec_pi     = np.array([0.5, 0.5])
    mean_rate  = len(spike_times) / (max_value - min_value)

    vec_lambda = np.empty(2)

    vec_lambda[0] = (mean_rate * 0.75) * bin_width
    vec_lambda[1] = (mean_rate * 1.25) * bin_width

    vec_spkt = np.array([spike - min_value for spike in spike_times])

    vec_Xi = get_vec_Xi(vec_spkt, bin_width)

    #########################################################
    #
    # バウムウェルチでのモデルパラメーターの最適化
    #
    ########################################################

    mat_Gamma, mat_Xi = HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)

    '''mat_A_old      = mat_A.copy()
    vec_pi_old     = vec_pi.copy()
    vec_lambda_old = vec_lambda.copy()
    なぜcopy()を使うのですか。代入より一桁多く時間がかかります'''
    mat_A_old = mat_A
    vec_pi_old = vec_pi
    vec_lambda_old = vec_lambda

    #####################
    # whileの評価について
    # モデルパラメータの変化が小さくなったところを終わりとしています
    # 具体的にはモデルパラメータの変化の差の和 sumcheck が
    # ある値より小さくなったときにflag=1として、
    # ループを抜けるようにしています
    # また、ループがあまりに多いときは抜けます
    #####################

    loop = 0
    flag = 0
    while(loop <= EMloop_num and flag == 0) :
        vec_pi_new, vec_lambda_new, mat_A_new = HMM_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi)

        '''vec_pi     = vec_pi_new.copy()
        vec_lambda = vec_lambda_new.copy()
        mat_A      = mat_A_new.copy()'''
        vec_pi     = vec_pi_new
        vec_lambda = vec_lambda_new
        mat_A      = mat_A_new

        sum_check = 0.0
        num_state = len(vec_pi)

        sum_check += sum(abs(vec_pi_old - vec_pi))

        sum_check += sum(abs(vec_lambda_old - vec_lambda))

        sum_check += sum(sum(abs(mat_A_old - mat_A)))

        if (sum_check / (1.0 * num_state * (num_state + 2)) < 1.0e-7) :
            flag = 1

        '''mat_A_old      = mat_A.copy()
        vec_pi_old     = vec_pi.copy()
        vec_lambda_old = vec_lambda.copy()'''
        mat_A_old      = mat_A
        vec_pi_old     = vec_pi
        vec_lambda_old = vec_lambda

        E_res     = HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)
        mat_Gamma = E_res[0]
        mat_Xi    = E_res[1]

        loop += 1

    print("loop",loop)
    #############################################
    #
    # ビタビアルゴリズムで最適な状態列を求める
    # stateは 0 or 1 で表される
    #
    #############################################

    vec_hidden = HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi)

    ############################################
    #
    # 以下で状態列の値をplotしやすいように,実際のrateに変換してる
    #
    ###########################################

    rate_func  = np.empty([len(vec_Xi), 2])

    c_time = 0.0

    for n in range(0, len(vec_Xi)) :
        state_id        = vec_hidden[n]
        rate_func[n][0] = round(c_time * 100) / 100.0
        rate_func[n][1] = round(vec_lambda[int(state_id)] * 100) / (bin_width * 100.0)

        c_time += bin_width

    return rate_func


def get_vec_Xi(vec_spkt, bin_width):
    """
    vec_spkt=初めのspikeがきた時を0としたspike時刻列
    bin_width=bin幅
    結果：
    観測列vec_Xiを出力。
    vec_Xiには各stepのspikeの個数が入ってます
    """
    spkt_dura = vec_spkt[len(vec_spkt) - 1]
    bin_num   = int(math.ceil(spkt_dura / bin_width))
    vec_Xi    = np.zeros(bin_num)

    ##########
    # spike時刻が各stepにあれば1を足していく
    ##########
    for x in vec_spkt:
        bin_id = int(math.floor(x / bin_width))
        if (bin_id < bin_num):
            vec_Xi[bin_id] += 1
    return vec_Xi


####
# HMM_E_step関数
# 隠れマルコフモデルの期待値を計算するステップです。

# 引数
# vec_Xi: ベクトルX_i、numpy arrayクラスで表現します。
# mat_A: 行列A、numpy arrayクラスで表現します。
# vec_lambda: ベクトルラムダ、numpy arrayクラスで表現します。
# vec_pi: ベクトルパイ、numpy arrayクラスで表現します。

# 返り値
# mat_Gamma:
# mat_Xi:
####

def HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi):
    """
    vec_Xi=観測列。各stepのspikeの個数が入ってます.
    mat_A=状態遷移確率行列
    vec_lambda=それぞれの状態でのbin内のspikeの平均値
    vec_pi=初期状態の確率
    """
    mat_emission      = get_mat_emission(vec_Xi, vec_lambda)
    vec_C, mat_alpha  = get_alpha_C(mat_A, vec_pi, mat_emission)
    mat_beta          = get_beta(mat_A, vec_pi, mat_emission, vec_C)
    mat_Gamma, mat_Xi = get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)

    '''res = [mat_Gamma, mat_Xi]

    return res'''
    return mat_Gamma, mat_Xi


def get_mat_emission(vec_Xi, vec_lambda):
    """
    vec_Xi=観測列
    vec_lambda=それぞれの状態でのbin内のspikeの平均値
    結果：
    stateごとに実際の各observationを得る確率を出力
    mat_emission=(step,state)の行列
    """
    ############
    # poisson分布を仮定しています。
    # vec_lambdaがpoisoon分布の平均値に対応します
    # stateごとに実際の各observationを得る確率を入れていきます
    #############
    mat_emission = np.array([[pow(x, y) * pow(math.e, -1.0 * x) / math.factorial(y) for x in vec_lambda] for y in vec_Xi])

    return mat_emission


def get_alpha_C(mat_A, vec_pi, mat_emission):
    """
    mat_A=状態遷移確率行列
    vec_pi=初期状態の確率
    mat_emission=時刻tでの観測を得る確率の行列.(step,state)
    結果：
    前向き変数alphaとスケーリングの係数vec_C(coefficient)を出力
    """

    num_of_states = len(vec_pi)
    num_of_obs = len(mat_emission)

    '''alpha_0 = np.array([mat_emission[0][i] * vec_pi[i] for i in range(0, num_of_states)])
    C_0     = 0.0
    for alpha in alpha_0 :
        C_0 += alpha

    vec_C_buf     = np.empty(num_of_obs)
    vec_C_buf[0]  = (C_0)
    alpha_0       = alpha_0 / C_0
    mat_alpha_buf = []
    mat_alpha_buf.append(list(alpha_0.copy()))

    for n in range(1, num_of_obs) :
        alpha_n = np.empty([num_of_states])
        for i in range(0, num_of_states) :
            sum_j = 0.0
            for j in range(0, num_of_states) :
                sum_j += mat_alpha_buf[n - 1][j] * mat_A[j][i]

            alpha_n_i = mat_emission[n][i] * sum_j
            alpha_n[i] = alpha_n_i

        C_n = 0.0
        for alpha in alpha_n :
            C_n += alpha

        vec_C_buf[n] = pd.Series([C_n])
        alpha_n = alpha_n / C_n

        mat_alpha_buf.append(list(alpha_n.copy()))

    res = [vec_C_buf, mat_alpha_buf]

    return res'''
    # 変更後
    alpha = np.empty((num_of_obs, num_of_states))
    coefficient = np.empty(num_of_obs)

    #####################
    # foward algorithm
    # initialization
    alpha[0, :] = vec_pi * mat_emission[0, :]
    # induction
    for t in range(0, num_of_obs - 1):  # scalingもしとく
        coefficient[t] = np.sum(alpha[t, :])
        alpha[t, :] = alpha[t, :] / coefficient[t]
        alpha[t+1, :] = (alpha[t, :] @ mat_A) * mat_emission[t+1, :]
    coefficient[num_of_obs - 1] = np.sum(alpha[num_of_obs - 1, :])
    alpha[num_of_obs - 1, :] = alpha[num_of_obs - 1, :] / coefficient[num_of_obs - 1]

    return coefficient, alpha


def get_beta(mat_A, vec_pi, mat_emission, vec_C):
    """
    mat_A=状態遷移確率行列
    vec_pi=初期状態の確率
    mat_emission=時刻tでの観測を得る確率の行列.(step,state)
    vec_C=alphaを計算するときに得られるスケーリング係数
    結果：
    後ろ向き変数betaを出力
    """
    num_of_states = len(vec_pi)
    num_of_obs    = len(mat_emission)

    # initialize1個目のデータCで割ってないやん
    mat_beta_buf = np.zeros([num_of_obs, num_of_states])

    for i in range(0, num_of_states) :
        mat_beta_buf[num_of_obs - 1][i] = 1.0

    for m in range(1, num_of_obs) :
        n               = num_of_obs - 1 - m
        mat_beta_buf_n1 = mat_beta_buf[n + 1]
        mat_emission_n1 = mat_emission[n + 1]
        mat_beta_buf_n  = mat_beta_buf[n]
        vec_C_n1        = vec_C[n + 1]
        for i in range(0, num_of_states) :
            sum_j = 0.0
            mat_A_i = mat_A[i]
            for j in range(0, num_of_states) :
                sum_j += mat_beta_buf_n1[j] * mat_emission_n1[j] * mat_A_i[j]

            mat_beta_buf_n[i] = (sum_j / vec_C_n1)

    return mat_beta_buf
    '''
    # 以下変更後 ここだけ元のプログラムの方が早かった。なぜだろう？
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
    mat_A=状態遷移確率行列
    mat_emission=時刻tでの観測を得る確率の行列.(step,state)
    mat_alpha=forward algorithmによって得られる変数
    mat_beta=backward algorithmによって得られる変数
    vec_C=alphaを計算するときに得られるスケーリング係数
    結果：
    gamma,xiを出力
    mat_Gamma=P(時刻tで状態i | vec_Xi, model)?
    mat_Xi=　?
    """
    num_of_states = len(mat_emission[0])
    num_of_obs    = len(mat_emission)

    '''mat_Gamma_buf = np.zeros([num_of_obs, num_of_states])いらない'''

    mat_Gamma_buf = mat_alpha * mat_beta

    '''mat_Xi_buf = np.zeros([num_of_obs - 1, num_of_states, num_of_states])
    for m in range(0, num_of_obs - 1) :
        mat_Xi_buf_m    = mat_Xi_buf[m]
        mat_alpha_m     = mat_alpha[m]
        mat_emission_m1 = mat_emission[m + 1]
        mat_beta_m1     = mat_beta[m + 1]
        vec_C_m1        = vec_C[m + 1]
        for i in range(0, num_of_states) :
            mat_Xi_buf_m_i = mat_Xi_buf_m[i]
            mat_alpha_m_i = mat_alpha_m[i]
            mat_A_i = mat_A[i]
            for j in range(0, num_of_states) :
                mat_Xi_buf_m_i[j] = (mat_alpha_m_i * mat_emission_m1[j] * mat_A_i[j] * mat_beta_m1[j]) / vec_C_m1
    なんでvec_Cでわるのかわかりません。
    betaを求める際にすでに割っていると思うのですが。
    下の変更はforを行列計算にしただけなので、アルゴリズムはかえてません'''

    xi = np.empty((num_of_obs - 1, num_of_states, num_of_states))
    for t in range(0, num_of_obs - 1):  # reshapeがあるので案外時間かかるかもしれない
        xi[t, :, :] = mat_alpha[t, :].reshape(
            num_of_states, 1) * mat_A * mat_emission[t + 1, :] * mat_beta[t + 1, :] / vec_C[t+1]

    '''res = [mat_Gamma_buf, mat_Xi_buf]'''

    return mat_Gamma_buf, xi

####
# HMM_M_step関数
# 隠れマルコフモデルの、尤度最大化ステップです。
####


def HMM_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi):
    """
    vec_Xi=観測列
    mat_A=状態遷移確率行列
    vec_lambda=それぞれの状態でのbin内のspikeの平均値
    vec_pi=初期状態の確率
    mat_Gamma=P(時刻tで状態i | vec_Xi, model)?
    mat_Xi= ?
    結果：
    一度更新したモデルパラメータ　vec_pi_new, vec_lambda_new, mat_A_new　を出力
    """
    num_of_states = len(mat_A)
    num_of_obs = len(vec_Xi)

    ####################
    # 新しいモデルパラメーターを求める
    ####################
    '''pi_denom = 0.0

    pi_denom += sum(mat_Gamma[0])

    vec_pi_new = mat_Gamma[0] / pi_denom
    まとめときます'''
    #################
    # pi
    ################
    vec_pi_new = mat_Gamma[0] / np.sum(mat_Gamma[0])

    # maxmize wrt lambda vector
    '''vec_lambda_new = np.empty(num_of_states)
    for k in range(0, num_of_states):
        lambda_denom = 0.0
        lambda_nume = 0.0
        for n in range(0, num_of_obs):
            lambda_denom += mat_Gamma[n][k]
            lambda_nume += mat_Gamma[n][k] * vec_Xi[n]

        if(lambda_denom == 0.0):
            vec_lambda_new[k] = 0.0
        else:
            vec_lambda_new[k] = lambda_nume / lambda_denom'''
    #################
    # vec_lambda
    #################
    vec_lambda_new = np.sum(mat_Gamma * vec_Xi.reshape(num_of_obs, 1), axis=0) / np.sum(mat_Gamma, axis=0)

    # maxmize wrt A matrix
    '''mat_A_new = np.zeros([num_of_states, num_of_states])
    for j in range(0, num_of_states) :
        A_denome = 0.0
        for n in range(0, num_of_obs - 1) :
            for mat_Xi_n_j_l in mat_Xi[n][j] :
                A_denome += mat_Xi_n_j_l

        for k in range(0, num_of_states) :
            A_nume = 0.0
            for mat_Xi_n in mat_Xi :
                A_nume += mat_Xi_n[j][k]
            if(A_denome == 0.0) :
                mat_A_new[j][k] = 0.0
            else :
                mat_A_new[j][k] = A_nume / A_denome'''
    ###############
    # a
    ###############
    total_gamma_a = np.zeros((1, num_of_states))
    for t in range(0, num_of_obs - 1):
        total_gamma_a += mat_Gamma[t, :]
    total_xi = np.sum(mat_Xi, axis=0)  # (2,2)の行列 = (i,j)
    mat_A_new = total_xi / total_gamma_a.T  # i同士の割り算

    '''res = [vec_pi_new, vec_lambda_new, mat_A_new]
    return res
    そのまま書いた方が代入ぶんの時間を省ける'''
    return vec_pi_new, vec_lambda_new, mat_A_new

def HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi):
    """
    vec_Xi=観測列
    mat_A=状態遷移確率行列
    vec_lambda=それぞれの状態でのbin内のspikeの平均値
    vec_pi=初期状態の確率
    結果：
    ビタビアルゴリズムで最適なstate列 vec_hs_seq を得て、
    出力する
    """
    mat_emission  = get_mat_emission(vec_Xi, vec_lambda)
    num_of_states = len(mat_A)
    num_of_obs    = len(vec_Xi)
    mat_hs_seq    = np.zeros([num_of_states, num_of_obs])
    vec_logp_seq  = np.zeros(num_of_states)

    for j in range(0, num_of_states) :
        mat_hs_seq[j][0] = j
        if (vec_pi[j] * mat_emission[0][j] == 0) :
            vec_logp_seq[j] = -np.inf
        else :
            vec_logp_seq[j]  = math.log(vec_pi[j] * mat_emission[0][j]) / math.log(10)

    for n in range(1, num_of_obs) :
        # copy the seq. up to n - 1
        mat_hs_seq_buf   = mat_hs_seq.copy()
        vec_logp_seq_buf = vec_logp_seq.copy()

        for j in range(0, num_of_states) :
            vec_h_logprob_i = np.zeros(num_of_states)
            for i in range(0, num_of_states) :
                vec_h_logprob_i[i] = vec_logp_seq[i] + math.log(mat_emission[n][j] * mat_A[i][j]) / math.log(10)

            max_element = max(vec_h_logprob_i)
            max_pos     = np.where(vec_h_logprob_i == max_element)[0][0]

            vec_logp_seq_buf[j] = max_element
            mat_hs_seq_buf[j]   = mat_hs_seq[max_pos].copy()

            mat_hs_seq_buf[j][n] = j

        mat_hs_seq   = mat_hs_seq_buf.copy()
        vec_logp_seq = vec_logp_seq_buf.copy()

    max_element = max(vec_logp_seq)

    max_pos = np.where(vec_logp_seq == max_element)[0][0]

    vec_hs_seq = mat_hs_seq[max_pos].copy()

    return vec_hs_seq


if __name__ == "__main__":
    # 実行するもの

    data = np.loadtxt("data4.txt")

    HMM(data)

