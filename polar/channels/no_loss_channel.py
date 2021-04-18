# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: no_loss_channel.py
Description: 极化码无损信道
======================================================
"""

from polar.channels.channel import Channel


class NoLossChannel(Channel):
    """无损信道类"""

    def __init__(self):
        pass

    def get_suffix(self):
        pass

    def get_erasure_prob(self):
        pass

    def get_llr(self, out_symbol):
        pass

    def get_ber(self):
        pass

    def modulate(self, to_message):
        pass

    def demodulate(self, from_message):
        pass

    def transmit(self, message):
        pass
