# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: bpsk_awgn_channel.py
Description: 实现BPSK + AWGN信道
======================================================
"""

import numpy as np
from scipy.special import erf

from polar.channels.channel import Channel

"""
===
相关参数：
signal noise ratio
SNR(dB)=10*log10(snr)

LR(Likelihood Ratio)

码率 R=K/N


_SINR_dB：以dB为单位的信噪比
_SINR：纯信噪比
_BER：比特误码率
_erasure_prob：删除概率 就是eposon
_zero_LLR：0的极大似然对数比
_one_LLR：1的极大似然对数比

这些参数是衡量一种解码方案乃至仿真系统性能的重要指标:
BER(bit error ratio): 比特误码率
FER(frame error ratio): 误帧率
BLER(block error ratio): 误块率
PER(package error ratio): 包错误率或丢包率
通过以SNR(dB)为横坐标，上述参数为纵坐标绘制折线图，
我们可以很直观的对比不同解码方案的性能

===

"""
class BpskAwgnChannel(Channel):
    """Bpsk_Awgn 信道"""
    def __init__(self, SINR, is_dB=True):
        """
        :description: Bpsk_Awgn信道构造函数
        :param 
            SINR(float): 信道比
            is_dB(bool): 信噪比单位是否为dB
            BER(float): 误码率
            zero_LLR(float): 0的似然对数比
            one_LLR(float)：1的似然对数比
        :Returns: 
        """
        super().__init__()
        if is_dB:
            self._SINR_dB = SINR
            self._SINR = np.power(10.0, SINR / 10.0)
        else:
            self._SINR = SINR
            self._SINR_dB = 10 * np.log10(SINR)
        # scipy.special使用scipy.special.erf()计算高斯曲线的面积
        # 误码率根据高斯参数求出来，使用密度函数
        self._BER = (1 - erf(np.sqrt(self._SINR))) / 2

        self._erasure_prob = - \
            (self._BER * np.log(self._BER) +
             (1.0 - self._BER) * np.log(1.0 - self._BER))

        # 0的似然值
        self._zero_LLR = np.log((1 - self._BER) / self._BER)
        # 1的似然值
        self._one_LLR = np.log(self._BER / (1 - self._BER))

    def get_suffix(self):
        """

        :return:
        """
        return 'SINR={}_dB'.format(self._SINR_dB)

    def get_erasure_prob(self):
        """

        :return:
        """
        return self._erasure_prob

    def get_llr(self, out_symbol):
        """

        :param out_symbol:
        :return:
        """
        return self._one_LLR if out_symbol == 1 else self._zero_LLR

    def get_ber(self):
        """

        :return:误码率
        """
        return self._BER

    def modulate(self, to_message):
        """
        调制
        :param to_message:
        :return:
        """
        return np.array([-1.0 if bit == 0 else 1.0 for bit in to_message], dtype='float64')

    def demodulate(self, from_message):
        """
        解调
        :param from_message:
        :return:
        """
        return np.array([0 if symbol <= 0.0 else 1 for symbol in from_message], dtype='uint8')

    def transmit(self, message):
        """
        传输
        :param message:
        :return:
        """
        noise_std = 1 / np.sqrt(2 * self._SINR)
        noise = noise_std * np.random.randn(len(message))

        return np.array(message + noise, dtype='float64')
    
    def __str__(self):
        s = """Bpsk_Awgn信道
        误块率: {}
        删除概率: {}
        One LLR: {}
        Zero LLR: {}
        suffix: {}
        """.format(self.get_ber(),self.get_erasure_prob(),
        self.get_llr(1),self.get_llr(0),self.get_suffix())
        return s


