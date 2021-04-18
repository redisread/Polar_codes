# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: matrix.py
Description: 极化码相关矩阵
======================================================
"""

import numpy as np
import math


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



if __name__ == "__main__":
    res = PolarMatrix.get_GN(4)
    print(res)
    # a = np.array([
    #     [1, 2],
    #     [3, 4]
    # ])
    # print(a[0,1])
    # print(PolarMatrix.get_GN(4))