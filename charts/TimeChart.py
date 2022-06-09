import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

f = open("../data/time.txt", "r")

time = [float(x.strip()) for x in f]
m = [i for i in range(1, len(time) + 1)]
f.close()

f = open("../data/time_fft.txt", "r")
time_fft = [float(x.strip()) for x in f]
m_fft = [i for i in range(1, len(time_fft) + 1)]
f.close()

plt.plot(m, time, label="Todd")
plt.plot(m_fft, time_fft, label="ToddFFT")
plt.xlabel("Todd degree")
plt.ylabel("Time")
plt.title("C++. Time dependence on degree")
plt.legend()
plt.show()