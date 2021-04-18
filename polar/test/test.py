# -*- coding:utf-8 -*-
"""
======================================================
Author: wujiahong
Filename: test.py
Description: 测试数学库
======================================================
"""
from polar.libcommon.matrix import PolarMatrix
from polar.code.encode import polar_encode_GN
from polar.libcommon.compute import get_sequence
from polar.libcommon.transmit import *
import numpy as np
# import sys
# sys.path.append(r"F:\github\Polar-Codes\polar_code")
# 添加路径


class Test(object):
    """测试类"""
    @classmethod
    def test_kron(cls) -> None:
        a = np.array([[1, 0], [1, 1]])
        r = np.kron(a, a)
        print(r)

    @classmethod
    def test_RN(cls):
        a = np.arange(0, 4).reshape((4, 1))
        print(a)
        RN = PolarMatrix.get_RN(4)
        print(RN)
        res = RN.dot(a)
        print(res)

    @classmethod
    def test_GN(cls):
        G8 = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                       [1, 0, 0, 0, 1, 0, 0, 0],
                       [1, 0, 1, 0, 0, 0, 0, 0],
                       [1, 0, 1, 0, 1, 0, 1, 0],
                       [1, 1, 0, 0, 0, 0, 0, 0],
                       [1, 1, 0, 0, 1, 1, 0, 0],
                       [1, 1, 1, 1, 0, 0, 0, 0],
                       [1, 1, 1, 1, 1, 1, 1, 1]], dtype=np.float64)
        g8 = PolarMatrix.get_GN(8)
        if (G8 == g8).all():
            print("GN Correct!")
        else:
            print("GN Error!")


def demodulate(X_message):
    return np.array([0 if symbol <= 0.0 else 1 for symbol
                     in X_message], dtype='uint8')


if __name__ == "__main__":
    N = 8
    K = 4
    epson = 0.5
    # message = np.random.randint(0, 2, K)
    message = np.array([0, 1, 1, 0])
    # message = np.ones(K)
    print("用户 message:", message)
    X_message = polar_encode_GN(message, K, N, epson)
    print("encode_message:", X_message)
    X_bpsk = bpsk(X_message)
    print("Bpsk message: ", X_bpsk)
    awgn_message = awgn(X_bpsk, 5)
    print("Awgn message:", awgn_message)

    print("Message to decode: ", demodulate(awgn_message))
