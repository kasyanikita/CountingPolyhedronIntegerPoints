import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

f = open("../data/bernoulli_time.txt", "r")

time = [float(x.strip()) for x in f]
f.close()

plt.plot(time, label="Bernoulli")
plt.xlabel("Bernoulli number")
plt.ylabel("Time")
plt.show()