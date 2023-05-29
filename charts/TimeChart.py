import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

f = open("data/time_fixed_det.txt", "r")
time = [float(x.strip()) for x in f]
f.close()
plt.figure(figsize=(13, 7))
plt.xticks(range(19), range(2, 21))
plt.plot(time)
plt.xlabel("Dimension")
plt.ylabel("Time (sec.)")
plt.title("Time dependence on dimension")
# plt.savefig("charts/pictures/TimeFixedDet.png")
plt.clf()

f = open("data/time_fixed_dim.txt", "r")
time = [float(x.strip()) for x in f]
f.close()
plt.figure(figsize=(13, 7))
plt.xticks(range(19, 300, 20), list(range(20, 301, 20)))
plt.plot(time)
plt.xlabel("Determinant value")
plt.ylabel("Time (sec.)")
plt.title("Time dependence on determinant")
# plt.show()
plt.savefig("charts/pictures/TimeFixedDim.png")