# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: 
Description: 计算相关
======================================================
"""

import numpy as np
import math

from polar.libcommon.matrix import PolarMatrix


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


def get_L1(y, N, v):
    """
    :description: 
    :param 
        y(int): 接受码
        N(int): 标号
        v(float): 方差
    :Returns: 
    """


def L_f(r1, r2):
    # y = (a*b + 1)/(a + b);%LR
    # y = log( (1+exp(a+b)) / (exp(a)+exp(b)) );%LLR
    # 硬件友好型LLR
    y = np.sign(r1) * np.sign(r2) * min(abs(r1), abs(r2))
    return y


def L_g(r1, r2, b):
    #% y = (a^(1 - 2*c)) * b;%LR
    #y = ( (-1)^c )*a + b;%LLR
    #
    y = r2 + (1 - 2*b)*r1
    return y

def judge(x):
    if x >= 0:
        return 0
    return 1

import numpy as np

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


class PolarCode(object):
    """极化码类"""
    N = None
    K = None

    def __init__(self, N: int, K: int):
        self.N = N
        self.K = K

    def encode(self):
        pass


if __name__ == "__main__":

    N = 4
    K = 3

    # 发送消息
    message = np.ones(K)
    # 获取极化序列
    res = get_sequence(np.log2(N)+1, 0.5)

    print("极化序列: ", res)

    # exit(1)
    # 获取比特混合下标
    r = np.argsort(res)
    print(r)
    # 输出
    print(res[r])

    # 初始化
    polarcode = PolarCode(N, K)

    # 获取混合比特消息
    u = np.zeros(N)
    # print(u[r[0:K]])
    u[r[0:K]] = message
    print("u: ", u)

    GN = PolarMatrix.get_GN(N)
    X = u.dot(GN) % 2

    print("x: ", X)

    # np.argpartition()
