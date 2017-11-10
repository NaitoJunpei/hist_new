import matplotlib.pyplot as plt
import math

def SS(spike_times) :
    max_value = max(spike_times)
    min_value = min(spike_times)
    onset     = min_value - 0.001 * (max_value - min_value)
    offset    = max_value + 0.001 * (max_value - min_value)

    for bin_num in range(1, 500) :
        bin_width = (offset - onset) / bin_num
        count     = [0.0] * bin_num
        for x in spike_times:
            count[int(math.floor((x - onset) / bin_width))] += 1

        av = 0.0
        va = 0.0
        for x in count :
            av += x / bin_num
            va += pow(x, 2.0) / bin_num

        cost = (2.0 * av - (va - av * av)) / pow(bin_width, 2.0)
        if (bin_num == 1 or cost < cost_min):
            # print(bin_num)
            cost_min        = cost
            optimal_bin_num = bin_num

    print(optimal_bin_num)
    drawSS(spike_times, optimal_bin_num)
    return optimal_bin_num

def drawSS(spike_times, optimal_bin_num):
    plt.hist(spike_times, optimal_bin_num)
    plt.show()
