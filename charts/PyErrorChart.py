import matplotlib.pyplot as plt
import sys
from os.path import dirname
import random
import seaborn as sns
sns.set_theme()

sys.path.insert(0, dirname(dirname(__file__)))

from PyTodd.Todd import Todd, ToddFFT

def todd3(n):
    v = [random.randint(1, 100) for _ in range(n)]
    triple_sum = 0
    for i in range(0, n - 2):
        for j in range(i + 1, n - 1):
            for k in range(j + 1, n):
                triple_sum += v[i] * v[j] * v[k]

    quad_pair_sum = 0
    for i in range(0, n):
        for j in range(0, n):
            if (i != j):
                quad_pair_sum += v[i] * v[i] * v[j]

    return v, quad_pair_sum / 24 + triple_sum / 8


error = []
errorFFT = []
n = 100
for i in range(1, n):
    v, ans = todd3(i)
    t = Todd(3, v)
    tfft = ToddFFT(3, v)
    t.main_todd()
    tfft.main_todd()
    error.append(abs(t.todd[3] - ans))
    errorFFT.append(abs(tfft.todd[3] - ans))
plt.plot(error, label="Todd")
plt.plot(errorFFT, label="ToddFFT")
plt.title("Python. Difference between expected and obtained values.")
plt.xlabel("Input size")
plt.ylabel("Error")
plt.legend()
plt.show()