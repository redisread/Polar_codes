# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: bsc_channel.py
Description: 实现BSC信道
======================================================
"""

import numpy as np

from polar.channels.channel import Channel

class BscChannel(Channel):
    """BSC信道，非对称信道"""

    def __init__(self, BER):
        super().__init__()
        self._BER = BER

        self._erasure_prob = - (self._BER * np.log(self._BER) + (1.0 - self._BER) * np.log(1.0 - self._BER))

        self._zero_LLR = np.log((1 - self._BER) / self._BER)
        self._one_LLR = np.log(self._BER / (1 - self._BER))

    def get_suffix(self):
        """

        :return:
        """
        return 'BER={}'.format(self._BER)

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

        :return:
        """
        return self._BER

    def modulate(self, to_message):
        """

        :param to_message:
        :return:
        """
        return np.array(to_message, dtype='uint8')

    def demodulate(self, from_message):
        """

        :param from_message:
        :return:
        """
        return np.array(from_message, dtype='uint8')

    def transmit(self, message):
        """

        :param message:
        :return:
        """
        error_pattern = np.array([0 if np.random.random_sample() > self._BER else 1 for _ in range(0, len(message))])

        return np.array(message ^ error_pattern, dtype='uint8')
