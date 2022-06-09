import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

f = open("../data/error.txt", "r")

error = [float(x.strip()) for x in f]
f.close()

f = open("../data/error_fft.txt", "r")
error_fft = [float(x.strip()) for x in f]
f.close()

plt.plot(error, label="Todd")
plt.plot(error_fft, label="ToddFFT")
plt.title("C++. Difference between expected values and obtained ones.")
plt.xlabel("Input size")
plt.ylabel("Error")
plt.legend()
plt.show()