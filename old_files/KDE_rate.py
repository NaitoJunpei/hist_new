##########
# KDE_rate.pyを実行するには、matplotlib、numpy、pandasライブラリが必要です

# 使い方
# KDE_rate.pyを、パスが通っているフォルダに置き、
# import KDE_rate
# をすると、ファイル内の関数が、KDE_rate.(関数名)の形で実行可能になります。

# ユーザーが使用するのはKDE関数のみで十分です。
# KDE関数は、spike列を引数に取ります。
# spike列の形式は、list、numpy.arrayなどが利用可能です。
# カーネル関数が最小となるビン幅を選び、カーネル関数をかけたグラフを描画します。
# 値を返しません。
##########

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

def KDE(spike_times) :
    spike_times = pd.Series(spike_times)

    optw = search_minimum(spike_times)
    opty = kernel(spike_times, optw)

    drawKDE(opty)

########## 
# search_minimum関数
# コスト関数が最小となるwを探し、その値を返します

# 引数
# spike_times: スパイク列

# 返り値
# コストが最小となるwの値
########## 

def search_minimum(spike_times) :
    max_value   = max(spike_times)
    min_value   = min(spike_times)
    onset       = min_value - 0.001 * (max_value - min_value)
    offset      = max_value - 0.001 * (max_value - min_value)
    
    def Cost(spike_times, w) :
        A = 0
        for i in range(0, len(spike_times)) :
            var1 = spike_times[i]
            for var2 in spike_times[i + 1:] :
                x = var1 - var2
                if (x < 5 * w) :
                    A += 2 * pow(math.e, (-x * x / (4 * w * w))) - 4 * math.sqrt(2) * pow(math.e, (-x * x / (2 * w * w)))

        return ((len(spike_times) + A) / (w * 2 * math.sqrt(math.pi)))
    
    C_min = np.inf
    for i in range(0, 50):
        W = (max_value - min_value) / (i + 1)
        C = Cost(spike_times, W)

        if (C < C_min) :
            C_min = C
            w = W

    return w

########## 
# kernel関数
# スパイク列に、設定したビン幅のGaussian Kernelをかけて返します

# 引数
# spike_times: スパイク列
# w: ビン幅

# 返り値
# spike_timesにGaussian Kernelをかけた値
##########

def kernel(spike_times, w) :
    max_value = max(spike_times)
    min_value = min(spike_times)
    K         = 200
    x         = xaxis(K, max_value, min_value)
    y         = [0] * K

    for i in range(0, K) :
        temp = 0
        for spike_time in spike_times :
            diff = x[i] - spike_time
            if (abs(diff) < 5 * w) :
                temp += gauss(diff, w)

        y[i] = temp

    return y

########## 
# gauss関数
# 平均0、分散w * wのgauss関数を計算します。
##########

def gauss(x, w) :
    return 1 / (w * math.sqrt(2 * math.pi)) * pow(math.e, (-x * x / (2 * w * w)))

##########
# xaxis関数
# グラフを描画する際のx方向成分の決定をします。

# 引数
# K: グラフを書く点の数
# max_value: x方向の最大値
# min_value: x方向の最小値

# 返り値
# x方向成分の配列
########## 

def xaxis(K, max_value, min_value) :
    x    = [0] * K
    x[0] = min_value
    for i in range(1, K) :
        x[i] = x[i - 1] + (max_value - min_value) / (K - 1)

    return x

########## 
# drawKDE関数
# 点列を受け取り、折れ線グラフの描画を行います。

# 引数
# opty: 描画する点列

# 返り値
# なし
##########

def drawKDE(opty) :
    plt.stackplot(range(0, 200), opty)
    plt.ylim(ymin = 0)
    plt.show()
