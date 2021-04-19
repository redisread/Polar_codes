# -*- coding:utf-8 -*-
"""
======================================================
Author: 吴嘉鸿(VictorHong)
Filename: polar_code.py
Description: 实现实现极化码编码译码类封装
======================================================
"""

import numpy as np
import matplotlib.pyplot as plt

from polar.code.construct import Construct
from polar.libcommon.utils import PolarMatrix
from polar.libcommon.utils import get_sequence
from polar.channels.bec_channel import BecChannel
from polar.channels.bpsk_awgn_channel import BpskAwgnChannel

# 算法4：遗传算法
# INPUT：N：码长，SNR：信道比
# OUTPUT：码长N的极化序列
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

        # 初始化极化码类
        self.snr = SNR
        self.channel = BpskAwgnChannel(self.snr, False)
        self.pc = PolarCode(self.n, self.K, self.channel, "PW")

        # 随机生成种群极化序列
        self.pop = np.vstack([np.random.permutation(self.DNA_size)
                              for _ in range(self.pop_size)])

    def translateDNA(self, DNA):
        """翻译DNA，直接返回极化下标序列"""
        return DNA

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
            A = self.pc.do_construct(self.pop[child])
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
            # select another individual from pop
            i_ = np.random.randint(0, self.pop_size, size=1)
            cross_points = np.random.randint(0, 2, self.DNA_size).astype(
                np.bool)   # choose crossover points
            # find the city number
            keep_DNA = parent[~cross_points]
            swap_DNA = pop[i_, np.isin(pop[i_].ravel(), keep_DNA, invert=True)]
            parent[:] = np.concatenate((keep_DNA, swap_DNA))
        return parent

    def mutate(self, child):
        for point in range(self.DNA_size):
            if np.random.rand() < self.mutate_rate:
                swap_point = np.random.randint(0, self.DNA_size)
                swapA, swapB = child[point], child[swap_point]
                child[point], child[swap_point] = swapB, swapA
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

    def compute(self, DNA):
        message = np.random.randint(0, 2, self.K)
        # message = np.ones(self.K)
        A = self.pc.do_construct(DNA)
        u, u_message = self.encode_decode(message)
        uc = (message == u_message[A])
        return uc.sum() / len(uc)

    def do_compute(self, DNA):
        message = np.random.randint(0, 2, self.K)
        # message = np.ones(self.K)
        A = self.pc.do_construct(DNA)
        u, u_message = self.encode_decode(message)
        uc = (message == u_message[A])

        # uc = (message == u_message[-self.K:])
        print("message: ", message)
        print("u: ", u)
        print("u_message: ", u_message)
        print(uc)
        return uc.sum() / len(uc)


def Genetic_Algorithm(N,SNR):
    """遗传算法"""
    POP_SIZE,CROSS_RATE,MUTATION_RATE= 50,0.2,0.02
    n = int(np.log2(N))
    R = 0.5
    K = int(N * R)
    GENERATIONS = 10
    Valid_t = 3
    pga = PolarGA(N, K, SNR, POP_SIZE, CROSS_RATE, MUTATION_RATE)
    # 开始迭代
    for generation in range(GENERATIONS):
        fitness = pga.get_fitness()
        best_DNA = pga.pop[np.argmax(fitness)]
        # print("Gen ", generation, " best DNA: ", best_DNA)
        rate = 0
        for _ in range(Valid_t):
            rate = rate + pga.compute(best_DNA)
        rate = rate / Valid_t
        # print("Rate:", rate)
        if rate == 1:
            break
        pga.evolve()
    return best_DNA


class PolarCode(object):
    """极化码类"""

    CRC_polynomials = {
        8: np.asarray([1, 1, 1, 0, 1, 0, 1, 0, 1], dtype='uint8'),
        16: np.asarray([1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], dtype='uint8'),
    }

    def __init__(self, n, K, channel, construction_method,CRC_len = 0):
        """
        :description: 构造函数
        :param  
                n(int): 定义码字长度，最终的长度为：(N = 2 ** n)
                N(int): 码字长度
                K(int): 信息比特的位数
                construction_methods(str): 极化码构造方法
                channel(Channel): 信道类型
        :Returns: 
        """

        self._n = n
        self._N = 2 ** n
        self._K = K
        self._channel = channel
        self._construction_method = construction_method
        self._K_minus_CRC = K
        self._CRC_len = CRC_len

        self._construction_methods = {
            'BathZ': Construct.Bhattacharyya,
            'MonteCalo': Construct.Monte_Carlo,
            'GA': Construct.GA,
            "Genetic": Genetic_Algorithm
        }

        # 解码方法字典
        self._decoding_methods = {
            'SC': self.sc_decode,
            'SCL': self._scl_decode,
        }

    def construct(self):
        """构造序列"""
        sequence = get_sequence(self._n+1, self._channel.get_erasure_prob())
        if self._construction_method not in self._construction_methods:
            sequence = Construct.Bhattacharyya(self._N,self._channel._SINR)
        elif self._construction_method == "MonteCalo":
            sequence = Construct.Monte_Carlo(self._N,self._channel._SINR,1000)
        else:
            sequence = self._construction_methods[self._construction_method](self._N,self._channel._SINR)
        sorted_sequence = np.argsort(sequence)
        A = sorted(sorted_sequence[0:self._K])
        A_c = sorted(sorted_sequence[self._K:])
        # A = sorted_sequence[0:self._K]
        # A_c = sorted_sequence[self._K:]
        # print("A(消息比特):", A, " A_c(冻结比特):", A_c)
        self._info_bits_positions = A
        self._frozen_bits_positions = A_c
        return A

    def do_construct(self, sequence):
        A = sorted(sequence[0:self._K])
        A_c = sorted(sequence[self._K:])
        # A = sequence[0:self._K]
        # A_c = sequence[self._K:]
        self._info_bits_positions = A
        self._frozen_bits_positions = A_c
        return A

    def fill_messgae(self, message):
        u = np.zeros(self._N)
        u[self._info_bits_positions] = message
        return u

    def encode(self, u):
        # u = np.zeros(self._N)
        # u[self._info_bits_positions] = message
        # print("混合信息比特: ", u)
        GN = PolarMatrix.get_GN(self._N)
        X = u.dot(GN) % 2
        return X

    def decode(self, y_message,way = 0):
        return self.sc_decode(y_message)
        # return self._scl_decode(y_message)

    def sc_decode(self, y_message):
        """连续消除译码"""
        y_message = np.array(y_message, dtype='uint8')
        u_est = np.zeros(self._N, dtype='uint8')
        for idx in self._info_bits_positions:
            llr = self.slow_llr(idx, self._N, y_message, u_est[:idx])
            u_est[idx] = 0 if llr > 0 else 1
        return u_est

    def bp_decode(self):
        """置信传播译码"""
        pass

    def ml_decode(self):
        """最大似然译码"""
        pass

    def slow_llr(self, i, N, y, u_est):
        # Trivial case of one polarized channel out of one.
        # 递归信道，直到N = 1才直接获取LLR值
        if i == 0 and N == 1:
            llr = self._channel.get_llr(y[i])

        else:
            if i % 2 == 0:
                # 上信道
                llr_1 = self.slow_llr(i // 2,
                                      N // 2,
                                      y[:(N // 2)],
                                      (u_est[::2] ^ u_est[1::2])[:(i // 2)])
                llr_2 = self.slow_llr(i // 2,
                                      N // 2,
                                      y[N // 2:],
                                      u_est[1::2][:(i // 2)])

                llr = llr_check_node_operation(llr_1, llr_2)
            else:
                # 下信道
                llr_1 = self.slow_llr((i - 1) // 2,
                                      N // 2,
                                      y[:(N // 2)],
                                      (u_est[:-1:2] ^ u_est[1:-1:2])[:(i - 1 // 2)])
                llr_2 = self.slow_llr((i - 1) // 2,
                                      N // 2,
                                      y[N // 2:],
                                      u_est[1::2][:((i - 1) // 2)])

                llr = llr_2 + ((-1) ** u_est[-1]) * llr_1

        return np.float64(llr)
    

    def _calculate_CRC(self, info_bits):
        """

        :param info_bits:
        :return:
        """
        if self._CRC_len == 0:
            return np.asarray([])
        else:
            padded_info_bits = np.concatenate([info_bits, np.zeros(self._CRC_len, dtype='uint8')])

            while len(padded_info_bits[0:self._K_minus_CRC].nonzero()[0]):
                cur_shift = (padded_info_bits != 0).argmax(axis=0)
                padded_info_bits[cur_shift: cur_shift + self._CRC_len + 1] ^= PolarCode.CRC_polynomials[self._CRC_len]

            return padded_info_bits[self._K_minus_CRC:]
    
    def _scl_decode(self, y_message, frozen_bits=None, L=32):
        self._initialize_data_structures(L)
        l = self._assign_initial_path()
        p_zero = self._get_array_pointer_p(0, l)

        for br in range(self._N):
            p_zero[br, 0] = self._out_bit_prob(y_message[br], 0)
            p_zero[br, 1] = self._out_bit_prob(y_message[br], 1)

        for phi in range(self._N):
            self._recursively_calc_p(self._n, phi, L)

            if phi in self._frozen_bits_positions:
                for l in range(L):
                    if self._active_path[l]:
                        c_curr = self._get_array_pointer_c(self._n, l)
                        c_curr[0, phi % 2] = frozen_bits[self._frozen_bits_positions.index(phi)] if frozen_bits is not None else 0
            else:
                self._continue_paths_unfrozen_bit(phi, L)

            if (phi % 2) == 1:
                self._recursively_update_c(self._n, phi, L)

        l_dash = 0
        p_dash = 0
        decoding_list = []
        is_CRC_present = False

        for l in range(L):
            if self._active_path[l]:
                path_output = self.polar_transform(self._get_array_pointer_c(0, l)[:, 0])
                path_output_info_bits = path_output[list(self._info_bits_positions)]

                if np.array_equal(self._calculate_CRC(path_output_info_bits[:self._K_minus_CRC]),
                                  path_output_info_bits[self._K_minus_CRC:]):
                    is_CRC_present = True
                    c_curr = self._get_array_pointer_c(self._n, l)
                    p_curr = self._get_array_pointer_p(self._n, l)
                    decoding_list.append(path_output)
                    if p_dash < p_curr[0, c_curr[0, 1]]:
                        l_dash = l
                        p_dash = p_curr[0, c_curr[0, 1]]

        if not is_CRC_present:
            return None

        c_zero = self._get_array_pointer_c(0, l_dash)
        return self.polar_transform(c_zero[:, 0])

    def _initialize_data_structures(self, L):
        """

        :param L:
        :return:
        """
        self._inactive_path_indices = []
        self._active_path = None
        self._array_pointer_p = [[] for _ in range(self._n + 1)]
        self._array_pointer_c = [[] for _ in range(self._n + 1)]
        self._path_index_to_array_index = None
        self._inactive_array_indices = [[] for _ in range(self._n + 1)]
        self._array_reference_count = None

        self._path_index_to_array_index = np.zeros((self._n + 1, L), dtype=np.uint8)
        self._array_reference_count = np.zeros((self._n + 1, L), dtype=np.uint8)

        for lam in range(self._n + 1):
            for s in range(L):
                self._array_pointer_p[lam].append(np.full(((2 ** (self._n - lam)), 2), -1.0))
                self._array_pointer_c[lam].append(np.zeros(((2 ** (self._n - lam)), 2), dtype=np.uint8))

                self._inactive_array_indices[lam].append(s)

        self._active_path = np.zeros(L, dtype=bool)
        for l in range(L):
            self._inactive_path_indices.append(l)
    
    def _assign_initial_path(self):
        """

        :return: Integer index of initial path
        """
        l = self._inactive_path_indices.pop()
        self._active_path[l] = True

        for lam in range(self._n + 1):
            s = self._inactive_array_indices[lam].pop()
            self._path_index_to_array_index[lam, l] = s
            self._array_reference_count[lam, s] = 1

        return l
    
    def _get_array_pointer_p(self, lam, l):
        """

        :param lam: An integer number of layer
        :param l: An integer path index
        :return: Reference to the corresponding probability pair array
        """

        s = self._path_index_to_array_index[lam, l]
        if self._array_reference_count[lam, s] == 1:
            s_dash = s
        else:
            s_dash = self._inactive_array_indices[lam].pop()
            self._array_pointer_p[lam][s_dash] = np.copy(self._array_pointer_p[lam][s])
            self._array_pointer_c[lam][s_dash] = np.copy(self._array_pointer_c[lam][s])
            self._array_reference_count[lam, s] -= 1
            self._array_reference_count[lam, s_dash] = 1
            self._path_index_to_array_index[lam, l] = s_dash

        return self._array_pointer_p[lam][s_dash]
    
    def _get_array_pointer_c(self, lam, l):
        """

        :param lam: An integer number of layer
        :param l: An integer path index
        :return: Reference to the corresponding bit pair array
        """

        s = self._path_index_to_array_index[lam, l]
        if self._array_reference_count[lam, s] == 1:
            s_dash = s
        else:
            s_dash = self._inactive_array_indices[lam].pop()
            self._array_pointer_c[lam][s_dash] = np.copy(self._array_pointer_c[lam][s])
            self._array_pointer_p[lam][s_dash] = np.copy(self._array_pointer_p[lam][s])
            self._array_reference_count[lam, s] -= 1
            self._array_reference_count[lam, s_dash] = 1
            self._path_index_to_array_index[lam, l] = s_dash

        return self._array_pointer_c[lam][s_dash]

    def _recursively_calc_p(self, lam, phi, L):
        """

        :param lam: An integer index of current layer
        :param phi: An integer index of current phase
        :param L: An integer size of decoding list
        :return: Void
        """

        if lam == 0:
            return
        psi = phi // 2

        if (phi % 2) == 0:
            self._recursively_calc_p(lam - 1, psi, L)

        sgm = 0.0
        for l in range(L):
            if self._active_path[l]:
                p_curr = self._get_array_pointer_p(lam, l)
                p_prev = self._get_array_pointer_p(lam - 1, l)
                c_curr = self._get_array_pointer_c(lam, l)

                for br in range(2 ** (self._n - lam)):
                    if (phi % 2) == 0:
                        p_curr[br, 0] = 0.5 * p_prev[2 * br, 0] * p_prev[2 * br + 1, 0] \
                                        + 0.5 * p_prev[2 * br, 1] * p_prev[2 * br + 1, 1]
                        sgm = max(sgm, p_curr[br, 0])

                        p_curr[br, 1] = 0.5 * p_prev[2 * br, 1] * p_prev[2 * br + 1, 0] \
                                        + 0.5 * p_prev[2 * br, 0] * p_prev[2 * br + 1, 1]
                        sgm = max(sgm, p_curr[br, 1])
                    else:
                        u = c_curr[br, 0]
                        p_curr[br, 0] = 0.5 * p_prev[2 * br, u] * p_prev[2 * br + 1, 0]
                        sgm = max(sgm, p_curr[br, 0])

                        p_curr[br, 1] = 0.5 * p_prev[2 * br, u ^ 1] * p_prev[2 * br + 1, 1]
                        sgm = max(sgm, p_curr[br, 1])

        for l in range(L):
            if self._active_path[l]:
                p_curr = self._get_array_pointer_p(lam, l)
                for br in range(2 ** (self._n - lam)):
                    p_curr[br, 0] /= sgm
                    p_curr[br, 1] /= sgm

    def _recursively_update_c(self, lam, phi, L):
        """

        :param lam: An integer index of current layer
        :param phi: An integer index of current phase
        :param L: An integer size of decoding list
        :return: Void
        """

        if (phi % 2) == 1:
            psi = phi // 2

            for l in range(L):
                if self._active_path[l]:
                    c_curr = self._get_array_pointer_c(lam, l)
                    c_prev = self._get_array_pointer_c(lam - 1, l)

                    for br in range(2 ** (self._n - lam)):
                        c_prev[2 * br][psi % 2] = c_curr[br][0] ^ c_curr[br][1]
                        c_prev[2 * br + 1][psi % 2] = c_curr[br][1]

            if psi % 2 == 1:
                self._recursively_update_c(lam - 1, psi, L)
    
    def _out_bit_prob(self, output_bit, input_bit):
        return self._channel.get_ber() if output_bit ^ input_bit else (1.0 - self._channel.get_ber())
    
    def _continue_paths_unfrozen_bit(self, phi, L):
        """

        :param phi: An integer index of current phase
        :param L: An integer size of decoding list
        :return: Void
        """
        prob_forks = np.zeros(2 * L, dtype=float)
        i = 0
        for l in range(L):
            if self._active_path[l]:
                p_curr = self._get_array_pointer_p(self._n, l)
                prob_forks[l] = p_curr[0, 0]
                prob_forks[l + L] = p_curr[0, 1]
                i += 1
            else:
                prob_forks[l] = -1
                prob_forks[l + L] = -1

        sorted_prob_forks = sorted(enumerate(prob_forks), key=lambda tup: -tup[1])
        rho = min(2 * i, L)
        cont_forks = np.zeros((L, 2), dtype=bool)
        for i in range(rho):
            cont_forks[sorted_prob_forks[i][0] % L, sorted_prob_forks[i][0] // L] = True
        for l in range(L):
            if self._active_path[l]:
                if not cont_forks[l][0] and not cont_forks[l][1]:
                    self._kill_path(l)

        for l in range(L):
            if not cont_forks[l][0] and not cont_forks[l][1]:
                continue
            c_curr = self._get_array_pointer_c(self._n, l)
            if cont_forks[l][0] and cont_forks[l][1]:
                c_curr[0][phi % 2] = 0
                l_dash = self._clone_path(l)
                c_curr = self._get_array_pointer_c(self._n, l_dash)
                c_curr[0][phi % 2] = 1
            elif cont_forks[l][0]:
                c_curr[0][phi % 2] = 0
            else:
                c_curr[0][phi % 2] = 1
    
    def _clone_path(self, l):
        """

        :param l: An integer index of path to clone
        :return: Integer index of copy
        """
        l_dash = self._inactive_path_indices.pop()
        self._active_path[l_dash] = True

        for lam in range(self._n + 1):
            s = self._path_index_to_array_index[lam, l]
            self._path_index_to_array_index[lam, l_dash] = s
            self._array_reference_count[lam, s] += 1

        return l_dash

    def _kill_path(self, l):
        """

        :param l: An integer index of path to kill
        :return:
        """
        self._active_path[l] = False
        self._inactive_path_indices.append(l)
        for lam in range(self._n + 1):
            s = self._path_index_to_array_index[lam, l]
            self._array_reference_count[lam, s] -= 1

            if self._array_reference_count[lam, s] == 0:
                self._inactive_array_indices[lam].append(s)
    
    @staticmethod
    def polar_transform(u_message):
        """
        Implements the polar transform on the given message in a recursive way (defined in Arikan's paper).

        :param u_message: An integer array of N bits which are to be transformed;
        :return: x_message -- result of the polar transform.
        """
        u_message = np.array(u_message)

        if len(u_message) == 1:
            x_message = u_message
        else:
            u1u2 = u_message[::2] ^ u_message[1::2]
            u2 = u_message[1::2]

            x_message = np.concatenate([PolarCode.polar_transform(u1u2), PolarCode.polar_transform(u2)])
        return x_message



def SGA(x):
    if x >= 0 and x < 0.036:
        y = 0.00206*x - 0.36*(x**2) + 18.68*(x**3)
    else:
        if x >= 0.036 and x < 0.05:
            y = -0.001288 + 0.05*x
        else:
            if x >= 0.05 and x < 0.92:
                y = 0.0159*x + 0.37*(x**2) - 0.112*(x**3)
            else:
                if x >= 0.92 and x < 8:
                    y = 0.367*x + 0.075*(x**2) - 0.0035*(x**3) - 0.18
                else:
                    y = -2.211 + 0.9848*x
    return y


def llr_check_node_operation(llr_1, llr_2):
    if llr_1 * llr_2 > 0:
        # If both LLRs are of one sign, we return the minimum of their absolute values.
        return min(abs(llr_1), abs(llr_2))
    else:
        # Otherwise, we return an opposite to the minimum of their absolute values.
        return -min(abs(llr_1), abs(llr_2))


class PolarTest(object):
    """极化码测试类"""
    @staticmethod
    def test_BEC():
        channel = BecChannel(0.5)
        pc = PolarCode(3, 4, channel, "PW")
        message = np.array([0, 1, 1, 0])
        pc.construct()
        u = pc.fill_messgae(message)
        X_message = pc.encode(u)
        print(X_message)

        Y_message = channel.transmit(X_message)
        print(Y_message)

    @staticmethod
    def encode_decode(pc,channel,u):
        X_message = pc.encode(u)

        X_message = channel.modulate(X_message)

        Y_message = channel.transmit(X_message)

        y_message = channel.demodulate(Y_message)

        u_message = pc.decode(y_message)

        return u_message
    
    @staticmethod
    def test_BPSK_AWGN(construct_method,N,K,SNR = 5):
        channel = BpskAwgnChannel(SNR, False)

        # BathZ MonteCalo GA Genetic
        pc = PolarCode(n, K, channel, construct_method)

        message = np.random.randint(0, 2, K)

        # 构造序列
        A = pc.construct()
        u = pc.fill_messgae(message)

        u_message = PolarTest.encode_decode(pc,channel,u)

        out_u = u_message[A]

        print((out_u == message).all())
    
    @staticmethod
    def test_BPSK_AWGN_log(construct_method,N,K,SNR = 5):
        channel = BpskAwgnChannel(SNR, False)

        # BathZ MonteCalo GA Genetic
        pc = PolarCode(n, K, channel, construct_method)

        message = np.random.randint(0, 2, K)

        print("消息内容：", message)

        # 构造序列
        A = pc.construct()

        u = pc.fill_messgae(message)
        print("填充比特：", u)

        X_message = pc.encode(u)
        print("编码信息比特：", X_message)

        X_message = channel.modulate(X_message)
        print("modulate:", X_message)

        Y_message = channel.transmit(X_message)
        print(Y_message)

        y_message = channel.demodulate(Y_message)
        print("dmodulate: ", y_message)

        u_message = pc.decode(y_message)
        print("Decode:", u_message)

        out_u = u_message[A]
        print("解码消息", out_u)
        print((out_u == message).all())
        print((message == u_message[-K:]).all())
        # out_u == message
    
    @staticmethod
    def test_BPSK_AWGN_BATCH(construct_method,N,K,T,SNR = 5):
        channel = BpskAwgnChannel(SNR, False)

        # BathZ MonteCalo GA Genetic
        n = int(np.log2(N))
        pc = PolarCode(n, K, channel, construct_method)

        # 构造序列
        A = pc.construct()
        c_sum = 0
        for _ in range(T):
            message = np.random.randint(0, 2, K)
            u = pc.fill_messgae(message)

            u_message = PolarTest.encode_decode(pc,channel,u)

            out_u = u_message[A]
            c_sum = c_sum + (out_u == message).all()
        return c_sum / T

def test_batch():
    SNR = 2
    channel = BpskAwgnChannel(1, False)
    n = 8
    N = 2 ** n
    R = 0.5
    K = int(N*R)
    # K = 3
    # BathZ MonteCalo GA
    pc = PolarCode(n, K, channel, "BathZ")
    M = 100
    sum = 0
    for _ in range(M):
        # message = np.array([1, 1, 1])
        message = np.random.randint(0, 2, K)
        # message = np.zeros(K)
        # print("消息内容：", message)
        # 构造序列
        A = pc.construct()
        # pc.do_construct([0,1,2,3])
        u = pc.fill_messgae(message)
        # print("填充比特：", u)

        X_message = pc.encode(u)
        # print("编码信息比特：", X_message)

        X_message = channel.modulate(X_message)
        # print("modulate:", X_message)

        Y_message = channel.transmit(X_message)
        # print(Y_message)

        y_message = channel.demodulate(Y_message)
        # print("dmodulate: ", y_message)

        u_message = pc.decode(y_message)
        # print("Decode:", u_message)

        out_u = u_message[A]
        # print("解码消息", out_u)
        # print((out_u == message).all())
        # correct = (out_u == message).all()
        correct = (u == u_message).all()
        sum = sum + correct
    print(sum/M)


def hhh():
    print("hhh")

if __name__ == "__main__":
    # n = 6
    # N = 2 ** n
    # R = 0.5
    # K = int(N*R)
    # SNR = 5
    # T = 100
    # # BathZ MonteCalo GA Genetic
    # # PolarTest.test_BPSK_AWGN("BathZ",N,K,SNR)

    # r = PolarTest.test_BPSK_AWGN_BATCH("Genetic",N,K,T,SNR)
    # print(r)


    cmethods = ["BathZ","MonteCalo","GA","Genetic"]
    test_N = 128
    test_R = 0.5
    test_K = int(test_N * test_R)
    test_T = 1

    M = 1000
    SNRs = [0.5 * i for i in range(1,10)]
    
    Rate = {}
    for method in cmethods:
        Rate[method] = np.empty(len(SNRs))
    
    y = np.empty_like(SNRs)
    for idx,snr in enumerate(SNRs):
        r = PolarTest.test_BPSK_AWGN_BATCH("Genetic",test_N,test_K,test_T,snr)
        y[idx] = 1 - r
    
    plt.plot(SNRs,y)
    plt.show()


    # for method in cmethods:
    #     idx = 0
    #     for snr in SNRs:
    #         r = PolarTest.test_BPSK_AWGN_BATCH(method,test_N,test_K,test_T,snr)
    #         Rate[method][idx] = 1 - r
    #         idx = idx + 1

    # for method in cmethods:
    #     plt.plot(SNRs,Rate[method],label=method)
    # plt.show()

    # test_batch()
    # test_BEC()
    # channel = BpskAwgnChannel(5, False)
    # pc = PolarCode(3, 4, channel, "PW")
    # seq = pc.GA()
    # print(seq)
