import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

f = open("../data/time_fixed_input.txt", "r")

time = [float(x.strip()) for x in f]
f.close()

f = open("../data/time_fft.txt", "r")
time_fft = [float(x.strip()) for x in f]
f.close()

plt.plot(time, label="Todd")
plt.plot(time_fft, label="ToddFFT")
plt.xlabel("Todd input size")
plt.ylabel("Time")
plt.title("Time dependence on input size")
plt.legend()
plt.show()