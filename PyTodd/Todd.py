from scipy.special import bernoulli
import matplotlib.pyplot as plt
import numpy as np
import time
from decimal import Decimal

class Todd:
    def __init__(self, _m, _xi):
        self.m = _m
        self.xi = _xi
        self.b = bernoulli(self.m)
        self.todd = [0] * (self.m + 1)
        self.todd_part = [0] * (self.m + 1)

    def update_part(self, x):
        self.todd_part[0] = 1
        fact = 1
        pow_x = 1
        for i in range(1, self.m + 1):
            fact *= i
            pow_x *= -x
            self.todd_part[i] = pow_x * self.b[i] / fact

    def init_todd(self):
        self.todd = self.todd_part[:]

    def calc_todd(self):
        res = [0] * (self.m + 1)
        for i in range(self.m + 1):
            for j in range(i + 1):
                res[i] += self.todd[i - j] * self.todd_part[j]
        return res

    def main_todd(self):
        self.update_part(self.xi[0])
        self.init_todd()
        for i in range(1, len(self.xi)):
            self.update_part(self.xi[i])
            self.todd = self.calc_todd()

class ToddFFT(Todd):
    def __init__(self, m, xi):
        super().__init__(m, xi)

    def calc_todd(self):
        n = 1
        while n < self.m + 1:
            n <<= 1
        n <<= 1
        self.todd += [0] * (n - self.m - 1)
        self.todd_part += [0] * (n - self.m - 1)
        fft_a = np.fft.rfft(self.todd)
        fft_b = np.fft.rfft(self.todd_part)
        for i in range(len(fft_a)):
            fft_a[i] *= fft_b[i]
        res = np.fft.irfft(fft_a)
        self.todd_part = self.todd_part[:self.m + 1]
        return list(res[:self.m + 1])
        