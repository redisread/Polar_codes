# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: draw.py
Description: 极化码极化过程演示
======================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import math


def phi(x):
    if x >= 0 and x <= 10:
        y = math.exp(-0.4527*x ^ 0.859 + 0.0218)
    else:
        y = math.sqrt(math.pi/x) * math.exp(-x/4) * (1 - 10/7/x)
    return y



def GA(N: int, sigma: float):
    n = np.log2(N)
    IWi = np.zeros((n+1, N))
    IWi[0, 0] = 1-epsilon
    for i in range(n):
        k = 2**i
        for j in range(k):
            tmp = IWi[i, j]
            IWi[i+1, 2*j-1] = tmp**2


def BEC(N, epsilon):
    index = int(np.log2(N))
    n = np.power(2, range(1, index+1))
    # 最大码长
    # max_code_length =
    W = np.zeros((n[index-1]+1, n[index-1]+1), dtype=np.float64)
    W[1][1] = epsilon

    for i in n:
        for j in range(1, int(i / 2)+1):
            W[i, 2 * j - 1] = W[int(i / 2), j] * W[int(i / 2), j]  # ^ 2;
            W[i, 2 * j] = 2 * W[int(i / 2), j] - \
                W[int(i / 2), j]*W[int(i / 2), j]   # ^ 2;

    # np.set_printoptions(threshold=np.inf)
    # print(W)
    plt.scatter(range(1, 2**index), W[2**index, 1:(2**index)], s=5, label='b.')
    plt.axis([0, n[index-1], 0, 1.0])
    plt.xlabel('Channel index')
    plt.ylabel('Symmetric capacity')
    plt.show()


if __name__ == '__main__':
    BEC(2**15, 0.5)
