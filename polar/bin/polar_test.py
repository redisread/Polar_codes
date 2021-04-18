# -*- coding:utf-8 -*-
"""
======================================================
Author: wujiahong
Filename: 
Description: polar code 测试入口
======================================================
"""

from polar.channels.bpsk_awgn_channel import BpskAwgnChannel
if __name__ == "__main__":
    channel = BpskAwgnChannel(5,is_dB=False)
    print(channel)

