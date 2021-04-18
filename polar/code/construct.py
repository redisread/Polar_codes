# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: construct.py
Description: 极化序列构造
======================================================
"""

import numpy as np
import random
from polar.libcommon.utils import PolarMatrix
from polar.libcommon.utils import get_sequence
from polar.channels.bec_channel import BecChannel
from polar.channels.bpsk_awgn_channel import BpskAwgnChannel
# from polar.code.polar_code import PolarCode

# 算法1：巴氏参数法
# INPUT：N：码长，SNR：信道比
# OUTPUT：码长N的极化序列


def bathZ(N, SNR):
    n = int(np.log2(N))
    z = np.zeros(N)
    z_init = np.exp(-SNR)
    z[0] = z_init
    for i in range(1, n+1):
        u = 2 ** i
        for t in range(0, u//2):
            T = z[t]
            z[t] = 2*T - T**2
            z[(u)//2 + t] = T**2
    res = np.zeros(N)
    for i in range(N):
        l = bit_reversed(i, n)
        res[l] = z[i]
    return res


# 相关API函数
def bit_reversed(x, n):
    """比特转换"""
    result = 0
    for i in range(n):  # for each bit number
        if (x & (1 << i)):  # if it matches that bit
            result |= (1 << (n - 1 - i))  # set the "opposite" bit in result
    return result


def UpdateL(n, L, i, j, B):
    """更新LLR"""
    u = 2**(n-j)
    l = i % u
    if l < u // 2:
        if L[i][j+1] == -1:
            UpdateL(n, L, i, j+1, B)
        if L[i+u//2][j+1] == -1:
            UpdateL(n, L, i+u//2, j+1, B)
        L[i][j] = (L[i][j+1]*L[i+u//2][j+1]+1)/(L[i][j+1]*L[i+u//2][j+1])
    else:
        if B[i-u//2][j] == 0:
            L[i][j] = L[i][j+1]*L[i-u//2][j+1]
        else:
            L[i][j] = L[i][j+1]/L[i-u//2][j+1]


def SGA(x):
    """近似高斯函数"""
    if x >= 0 and x < 0.036:
        y = 0.00206*x - 0.36*(x**2) + 18.68*(x**3)
    else:
        if x >= 0.036 and x < 0.05:
            y = -0.001288 + 0.05*x
        else:
            if x >= 0.05 and x < 0.92:
                y = 0.0159*x + 0.37*(x**2) - 0.112*(x**3)
            else:
                if x >= 0.92 and x < 8:
                    y = 0.367*x + 0.075*(x**2) - 0.0035*(x**3) - 0.18
                else:
                    y = -2.211 + 0.9848*x
    return y


# 算法2：蒙特卡洛法
# INPUT：N：码长，SNR：信道比，M：迭代次数
# OUTPUT：码长N的极化序列
def Monte_Carlo_estimation(N, SNR, M):
    """蒙特卡洛构造法"""
    n = int(np.log2(N))
    c = 0
    variance = 1 / SNR
    for i in range(M):
        # 传输全0码字
        y = np.sqrt(SNR) + np.random.randn(N)

        L = np.zeros((n+1, N))
        for j in range(N):
            L[0, j] = 2 * y[j]/variance

        for j in range(1, n+1):
            for k in range(0, N):
                if k % 2 == 0:
                    L[j, k] = SGA(L[j-1, k//2])
                else:
                    L[j, k] = 2 * L[j-1, (k-1)//2]
        L = L[n, :]
        d = np.zeros(N)
        for i in range(N):
            l = bit_reversed(i, n)
            d[l] = 0 if L[l] >= 1 else 1
        c = c + d
    return c/M


def phi(x):
    if x > 0 and x <= 7.0633:
        return np.exp(0.0116*(x**2) - 0.4212*x)
    else:
        return np.exp(-0.2944*x - 0.3169)


def phi_reverse(x):
    if x > 0 and x < 0.0911:
        return -(np.log(x) + 0.3169) / 0.2944
    else:
        return (0.4212 - np.sqrt(0.4212**2 + 4 * 0.0116*np.log(x))) / (2 * 0.0116)


def Gaussian_Approximation(N, SNR):
    """高斯近似构造法"""
    variance = 1 / SNR  # 计算方差
    n = int(np.log2(N))
    llr = np.zeros((n+1, N))
    llr[0, :] = 2 / variance
    for i in range(1, n + 1):
        for j in range(0, N):
            if j % 2 == 0:
                llr[i, j] = phi_reverse(1 - (1 - phi(llr[i-1, j // 2]))**2)
            else:
                llr[i, j] = 2 * llr[i-1, (j - 1)//2]
    return llr[n, :]

# 算法3：高斯近似法
# INPUT：N：码长，SNR：信道比
# OUTPUT：码长N的极化序列
def Gaussian_Approximation2(N, SNR):
    """高斯近似法计算LLR"""
    n = int(np.log2(N))
    variance = 1 / SNR
    llr = np.zeros((n+1, N))
    llr[0, :] = 2 / variance
    for i in range(1, n + 1):
        for j in range(0, N):
            if j % 2 == 0:
                llr[i, j] = SGA(llr[i-1, j // 2])
            else:
                llr[i, j] = 2 * llr[i-1, (j - 1)//2]
    return llr[n, :]



class Construct(object):
    """极化序列构造类"""
    # 构造一个k维的极化码即找到k个最可靠的位信道#
    @staticmethod
    def Bhattacharyya(N, SNR):
        """巴氏参数法"""
        return bathZ(N, SNR)

    @staticmethod
    def Monte_Carlo(N, SNR, M):
        """蒙特卡洛法"""
        return Monte_Carlo_estimation(N, SNR, M)

    @staticmethod
    def GA(N, SNR):
        """高斯近似法"""
        return Gaussian_Approximation2(N, SNR)[::-1]
    
    # @staticmethod
    # def Genetic(N,SNR):
    #     """遗传算法"""
    #     return Genetic_Algorithm(N,SNR)


def getSeq(seq):
    sorted_sequence = np.argsort(seq)
    return sorted_sequence


def test_bathZ(N, SNR):
    seq = Construct.Bhattacharyya(N, SNR)
    print(seq)
    print(getSeq(seq))


def test_GA(N, SNR):
    seq = Construct.GA(N, SNR)
    print(seq)
    print(getSeq(seq))


def test_MK(N, SNR, M):
    seq = Construct.Monte_Carlo(N, SNR, M)
    print(seq)
    print(getSeq(seq))


if __name__ == "__main__":
    # test
    n = 2
    N = 2 ** n
    M = 1000
    SNR = 5
    # test_bathZ(N,SNR)
    test_GA(N, SNR)
    # test_MK(N,SNR,M)
