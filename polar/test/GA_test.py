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

# import matplotlib.pyplot.plot

class PolarGA(object):
    """遗传算法极化码构造类"""

    def __init__(self, N, K, message,SNR, POP_SIZE, CROSS_RATE, MUTATE_RATE):
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

        self.message = message
        self.set_message = message
        self.sey_frozen = np.zeros(N-K)

        # 初始化极化码类
        self.snr = SNR
        self.channel = BpskAwgnChannel(self.snr, False)
        self.pc = PolarCode(self.n, self.K, self.channel, "PW")

        # 随机生成种群极化序列emvstack([np.random.permutationpop:DNA
        self.pop = np.zeros((self.pop_size,self.DNA_size),dtype=np.uint8)

        for i in range(self.pop_size):
            seq = np.random.permutation(self.DNA_size)
            to_seq = np.zeros(self.DNA_size)
            to_seq[seq[0:K]] = message
            self.pop[i,:] = to_seq 

    def translateDNA(self, DNA):
        """翻译DNA，直接返回极化下标序列"""
        return DNA

    def decode(self, DNA):
        """比特混合+编码解码"""
        y_message = DNA.copy()
        u_message = self.pc.decode(y_message)
        # 返回混合比特和解码比特
        return u_message

    def get_fitness(self):
        """获取适应度"""
        # 对每一个极化序列进行译码
        compare = np.zeros(self.pop_size)
        # 直接对所有的极化序列进行译码，不用构造
        for child in range(self.pop_size):
            u_message = self.decode(self.pop[child])
            uc = (self.message == u_message[:self.K])
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
                child[point] = np.random.randint(*self.DNA_bound)  # choose a random ASCII index
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


def test_PGA(N):
    n = int(np.log2(N))
    # N = 2 ** n
    SNR = 1
    R = 0.5
    K = int(N * R)
    message = np.random.randint(0,2,K)
    POP_SIZE = 20
    CROSS_RATE = 0.3
    MUTATE_RATE = 0.1

    pga = PolarGA(N,K,message,SNR,POP_SIZE,CROSS_RATE,MUTATE_RATE)

    generation = 10
    print("message:",message)
    for gen in range(generation):
        fitness = pga.get_fitness()
        best_DNA = pga.pop[np.argmax(fitness)]
        print("Gen ",gen," best DNA: ",best_DNA)
        if (best_DNA[:K] == message).all():
            break
        pga.evolve()
    print("decode:",best_DNA)



if __name__ == "__main__":
    n = 2
    N = 2 ** n
    test_PGA(N)
    # test_rate()
