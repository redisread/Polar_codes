# -*- coding:utf-8 -*-
"""
======================================================
Author: wujiahong
Filename: test.py
Description: 测试数学库
======================================================
"""

import numpy as np

if __name__ == "__main__":
    a = np.array([[1,0],[1,1]])
    r = np.kron(a,a)
    print(r)