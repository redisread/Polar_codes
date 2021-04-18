# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: bec_channel.py
Description: 实现BEC信道
======================================================
"""

import numpy as np
import random
from polar.channels.channel import Channel


class BecChannel(Channel):
    """BEC信道"""

    def __init__(self, SNR):
        """构造函数，参数为误码率，误码率等同于删除概率"""
        super().__init__()
        BER = 1 / (1 + SNR)
        self._BER = BER
        self._erasure_prob = BER
        self._SINR = SNR
        self._SINR_dB = 10 * np.log10(SNR)

    def get_erasure_prob(self):
        """

        :return:
        """
        return self._erasure_prob
    
    def transmit(self,message):
        """传输信息，随机删除一些比特，设置删除值为255"""
        return self.erase(message,self._erasure_prob)
    
    @staticmethod
    def erase(message, epson) -> np.ndarray:
        """随机删除操作"""
        out_message = np.array(message)
        length = len(message)
        sample_num = int(epson * length)
        pos_indexs = [i for i in range(length)]
        sample_list = random.sample(pos_indexs, sample_num)
        out_message[sample_list] = 255
        return out_message
    
    def get_llr(self, out_symbol):
        
        pass
        
    def get_ber(self):
        pass

    
    def modulate(self, to_message):
        pass

    
    def demodulate(self, from_message):
        pass
    def __str__(self):
        s = """BEC信道
        删除概率: {}
        One LLR: {}
        Zero LLR: {}
        suffix: {}
        """.format(self.get_erasure_prob(),
        self.get_llr(1),self.get_llr(0),self.get_suffix())
        return s
    
    def get_suffix(self):
        pass

if __name__ == "__main__":
    message = np.random.randint(0, 2, 16)
    m = BecChannel.erase(message,0.5)
    print(m)
