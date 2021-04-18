# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: 
Description: 计算相关
======================================================
"""

# from polar.libcommon.matrix import get_GN
from ..libcommon.compute import get_sequence
from ..libcommon.matrix import PolarMatrix
import numpy as np


def polar_encode(U, N, SNR):
    """
    :description: 极化码编码
    :param 
        U: 消息
        n: 码长
        SNR: 信噪比db
    :Returns: 
    """
    snr = 10**(SNR//10)

    # 挑选信息集合

    pass


def polar_encode_GN(message, K, N, epson):
    sequence = get_sequence(np.log2(N)+1, epson)
    sorted_sequence = np.argsort(sequence)
    A = sorted(sorted_sequence[0:K])
    A_c = sorted(sorted_sequence[K:])
    print("A(消息比特):", A, " A_c(冻结比特):", A_c)
    u = np.zeros(N)
    u[A] = message
    print("混合信息比特: ", u)
    GN = PolarMatrix.get_GN(N)
    X = u.dot(GN) % 2
    return X




if __name__ == "__main__":
    pass
