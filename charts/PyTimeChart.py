import matplotlib.pyplot as plt
import sys
from os.path import dirname
import time
import seaborn as sns
sns.set_theme()

sys.path.insert(0, dirname(dirname(__file__)))

from PyTodd.Todd import Todd, ToddFFT

times = []
timesFFT = []

for i in range(1, 100):
    t = Todd(i, [7, 10, 5, 4])
    start = time.time()
    t.main_todd()
    end = time.time()
    times.append(end - start)


for i in range(1, 100):
    t = ToddFFT(i, [7, 10, 5, 4])
    start = time.time()
    t.main_todd()
    end = time.time()
    timesFFT.append(end - start)

plt.plot(times, label="Todd")
plt.plot(timesFFT, label="ToddFFT")
plt.xlabel("Todd degree")
plt.ylabel("Time")
plt.title("Python. Time dependence on degree")
plt.legend()
plt.show()