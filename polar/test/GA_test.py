# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: GA_test.py
Description: 实现遗传算法测试
======================================================
"""

import numpy as np
import random
from polar.libcommon.utils import PolarMatrix
from polar.libcommon.utils import get_sequence
from polar.channels.bec_channel import BecChannel
from polar.channels.bpsk_awgn_channel import BpskAwgnChannel
from polar.code.polar_code import PolarCode
from polar.code.construct import Construct

# import matplotlib.pyplot.plot

class PolarGA(object):
    """遗传算法极化码构造类"""

    def __init__(self, N, K, SNR, POP_SIZE, CROSS_RATE, MUTATE_RATE):
        """
        :description: 初始化函数
        :param 
            N: 码长
            K: 信息比特数
            SNR: 信噪比
            POP_SIZE: 种群数
            CROSS_ROTE: 交叉概率
            MUTATE_RATE: 变异概率
        :Returns: 
        """

        self.n = int(np.log2(N))
        self.N = N
        self.K = K
        self.DNA_size = self.N
        self.pop_size = POP_SIZE

        self.cross_rate = CROSS_RATE
        self.mutate_rate = MUTATE_RATE


        self.NUM_SLICE = 1000
        self.SNR_MIN = 0.5
        self.SNR_MAX = 5
        self.SNR_SCLCE = (self.SNR_MAX - self.SNR_MIN) / self.NUM_SLICE
        self.DNA_size = int(np.log2(self.NUM_SLICE)) 


        # 初始化极化码类
        self.snr = SNR
        self.channel = BpskAwgnChannel(self.snr, False)
        self.pc = PolarCode(self.n, self.K, self.channel, "BathZ")

        # 随机生成种群极化序列
        self.pop = np.vstack([np.random.randint(0,2,self.DNA_size)
                              for _ in range(self.pop_size)])

    def get_SNR(self,DNA):
        n = len(DNA)
        num = 0
        for idx,val in enumerate(DNA):
            num = num + val * (2 ** (n - 1 - idx))
        return self.SNR_MIN + num * self.SNR_SCLCE


    def translateDNA(self, DNA):
        """翻译DNA，直接返回极化下标序列"""
        snr = self.get_SNR(DNA)
        seq = Construct.Bhattacharyya(self.N,snr)
        return np.argsort(seq)

    def encode_decode(self, message):
        """比特混合+编码解码"""
        # 混合比特
        u = self.pc.fill_messgae(message)
        X_message = self.pc.encode(u)
        X_message = self.channel.modulate(X_message)
        Y_message = self.channel.transmit(X_message)
        y_message = self.channel.demodulate(Y_message)

        # 解码比特
        u_message = self.pc.decode(y_message)
        # 返回混合比特和解码比特
        return u, u_message

    def get_fitness(self):
        """获取适应度"""
        # 对每一个极化序列进行译码
        # compare = np.zeros((self.pop_size,self.DNA_size))
        compare = np.zeros(self.pop_size)
        message = np.random.randint(0, 2, self.K)
        # 直接对所有的极化序列进行译码，不用构造
        for child in range(self.pop_size):
            seq = self.translateDNA(self.pop[child])
            A = self.pc.do_construct(seq)
            u, u_message = self.encode_decode(message)
            uc = (message == u_message[A])
            compare[child] = uc.sum()

        return compare / self.DNA_size

    def select(self, fitness):

        idx = np.random.choice(np.arange(self.pop_size),
                               size=self.pop_size,
                               replace=True,
                               p=fitness/fitness.sum())
        return self.pop[idx]

    def crossover(self, parent, pop):
        if np.random.rand() < self.cross_rate:
            i_ = np.random.randint(0, self.pop_size, size=1)                        # select another individual from pop
            cross_points = np.random.randint(0, 2, self.DNA_size).astype(np.bool)   # choose crossover points
            parent[cross_points] = pop[i_, cross_points]                            # mating and produce one child
        return parent

    def mutate(self, child):
        for point in range(self.DNA_size):
            if np.random.rand() < self.mutate_rate:
                child[point] = np.random.randint(0,self.DNA_size)  # choose a random ASCII index
        return child

    def evolve(self):
        fitness = self.get_fitness()
        pop = self.select(fitness)
        pop_copy = pop.copy()
        for parent in pop:  # for every parent
            child = self.crossover(parent, pop_copy)
            child = self.mutate(child)
            parent[:] = child
        self.pop = pop

    def compute(self, seq):
        message = np.random.randint(0, 2, self.K)
        A = self.pc.do_construct(seq)
        u, u_message = self.encode_decode(message)
        uc = (message == u_message[A])
        return uc.sum() / len(uc)
    
    def Valid_SNR(self,seq):
        sum = 0
        for _ in range(3):
            message = np.random.randint(0, 2, self.K)
            A = self.pc.do_construct(seq)
            u, u_message = self.encode_decode(message)
            sum =  sum + (message == u_message[A]).all()
        return sum / 3
    
    def do_compute(self, DNA):
        message = np.random.randint(0, 2, self.K)
        # message = np.ones(self.K)
        A = self.pc.do_construct(DNA)
        u, u_message = self.encode_decode(message)
        uc = (message == u_message[A])

        # uc = (message == u_message[-self.K:])
        print("message: ",message)
        print("u: ",u)
        print("u_message: ",u_message)
        print(uc)
        return uc.sum() / len(uc)
        


def test_PGA():
    POP_SIZE = 3
    CROSS_RATE = 0.1
    MUTATION_RATE = 0.02

    # 极化码长度的幂次
    n = 8
    # 极化码长度
    N = 2 ** n
    # 极化码的码率
    R = 0.5
    K = int(N * R)

    # 信噪比
    SNR = 1

    # 迭代的子代数量
    GENERATIONS = 20

    pga = PolarGA(N, K, SNR, POP_SIZE, CROSS_RATE, MUTATION_RATE)

    # 开始迭代
    for generation in range(GENERATIONS):
        fitness = pga.get_fitness()
        best_DNA = pga.pop[np.argmax(fitness)]
        design_snr = pga.get_SNR(best_DNA)
        best_seq = pga.translateDNA(best_DNA)
        print("Gen ", generation, "best design-snr：",
                design_snr," best Sequence: ", best_seq)
        rate = pga.Valid_SNR(best_seq)
        print("Rate:", rate)
        if rate == 1:
            break
        pga.evolve()

    print("验证误码率:")
    m = 50
    sum = 0
    for i in range(m):
        rate = pga.compute(best_seq)
        if rate == 1:
            sum = sum + 1
    print("Rate: ", sum / m)


def test_rate():
    POP_SIZE = 20
    CROSS_RATE = 0.1
    MUTATION_RATE = 0.02

    # 极化码长度的幂次
    n = 6
    # 极化码长度
    N = 2 ** n
    # 极化码的码率
    R = 0.5
    K = int(N * R)

    # 信噪比
    SNR = 5

    # 迭代的子代数量
    GENERATIONS = 20

    pga = PolarGA(N, K, SNR, POP_SIZE, CROSS_RATE, MUTATION_RATE)

    best_DNA = np.random.permutation(N)
    # best_DNA = [3,2,1,0]
    print("DNA: ",best_DNA)
    rate = pga.do_compute(best_DNA)
    print(rate)
    # if rate == 1:
    #     m = 50
    #     sum = 0
    #     for i in range(m):
    #         rate = pga.compute(best_DNA)
    #         if rate == 1:
    #             sum = sum + 1
    #     print("Rate: ", sum / m)


if __name__ == "__main__":
    test_PGA()
    # test_rate()
