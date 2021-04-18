# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: utils.py
Description: 计算相关
======================================================
"""

import numpy as np


def bpsk(X):
    """BPSK调制"""
    return 2*X - 1

def awgn(x, snr, seed=7):
    '''
    加入高斯白噪声 Additive White Gaussian Noise
    :param x: 原始信号
    :param snr: 信噪比
    :return: 加入噪声后的信号
    '''
    np.random.seed(seed)  # 设置随机种子
    snr = 10 ** (snr / 10.0)
    xpower = np.sum(x ** 2) / len(x)
    npower = xpower / snr
    noise = np.random.randn(len(x)) * np.sqrt(npower)
    return x + noise



def get_sequence(n: int, epson: float) -> np.ndarray:
    """
    :description: 极化
    :param n(int): 序列大小
    :Returns: 极化之后的序列
    考虑使用迭代的方式
    """

    if n == 1:
        sequence = np.array([epson])
    else:
        epson1 = 2 * epson - epson**2
        epson2 = epson**2
        sequence1 = get_sequence(n-1, epson1)
        sequence2 = get_sequence(n-1, epson2)
        sequence = np.hstack((sequence1, sequence2))
    return sequence


class PolarMatrix(object):
    """极化码矩阵类"""

    @classmethod
    def get_RN(cls, N: int) -> np.ndarray:
        """
        :description: 获取RN矩阵，实现输入S(1,n) 到 (s1,s3,...,s2,s4,...)
        :param N(int): 矩阵大小
        :Returns: RN
        """

        RN = np.zeros((N, N))
        for i in range(0, N):
            if i % 2 == 0:
                RN[i][i//2] = 1
            else:
                RN[i][(i+N-1)//2] = 1
        return RN

    @classmethod
    def get_BN(cls, N: int) -> np.ndarray:
        """
        :description: 获取BN矩阵，递归式 B_N = R_N(I_2 kron B_{N/2}) B2 = I2
        :param N(int): 矩阵的大小,必须的2的幂次
        :Returns: BN
        """

        if N == 2:
            BN = np.eye(2)
        else:
            RN = PolarMatrix.get_RN(N)
            BN = RN.dot(np.kron(np.eye(2), PolarMatrix.get_BN(N//2)))
        return BN

    @classmethod
    def get_GN(cls, N: int) -> np.ndarray:
        """
        :description: 获取GN矩阵,G_N = B_N * (N次F的kron积) 
        :param N(int): GN矩阵的大小
        :Returns: GN
        """
        n = int(np.log2(N))
        f = np.array([
            [1, 0],
            [1, 1]
        ])
        F = f
        for _ in range(0, n-1):
            F = np.kron(F, f)
        BN = PolarMatrix.get_BN(N)
        GN = BN.dot(F)
        return GN

def bit_reversed(x, n):
    result = 0
    for i in range(n):  # for each bit number
        if (x & (1 << i)):  # if it matches that bit
            result |= (1 << (n - 1 - i))  # set the "opposite" bit in result
    return result

if __name__ == "__main__":
    a = np.arange(4)
    BN = PolarMatrix.get_BN(4)
    print(a)
    print(BN)
    print(a.dot(BN))
