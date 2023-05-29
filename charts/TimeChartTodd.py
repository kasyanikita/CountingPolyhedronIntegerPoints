import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

f = open("data/time_fixed_degree.txt", "r")
time = [float(x.strip()) for x in f]
f.close()
plt.figure(figsize=(13, 7))
plt.plot(time, label="Todd")
plt.xlabel("Todd input size")
plt.ylabel("Time")
plt.title("Time dependence on input size")
plt.legend()
plt.savefig("charts/pictures/TimeFixedDegree.png")
plt.clf()

f = open("data/time_fixed_input.txt", "r")
time = [float(x.strip()) for x in f]
f.close()

f = open("data/time_fft_fixed_input.txt", "r")
time_fft = [float(x.strip()) for x in f]
f.close()

plt.plot(time, label="Todd")
plt.plot(time_fft, label="ToddFFT")
plt.xlabel("Todd degree")
plt.ylabel("Time")
plt.title("Time dependence on degree")
plt.legend()
plt.savefig("charts/pictures/TimeFixedInput.png")