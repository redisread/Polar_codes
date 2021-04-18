# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: 
Description: 
======================================================
"""

import numpy as np
from ..libcommon.compute import get_sequence


def SGA(x):
    if x >= 0 and x < 0.036:
        y = 0.00206*x - 0.36*x^2 + 18.68*x^3


def decode_GN(X,N,K,epson):
    # 挑选信息集
    sequence = get_sequence(np.log2(N)+1, epson)
    sorted_sequence = np.argsort(sequence)
    A = sorted_sequence[0:K]
    A_c = sorted_sequence[K:]
    tmp = np.zeros(N)
    # tmp[A] = message

    pass

if __name__ == "__main__":
    pass
