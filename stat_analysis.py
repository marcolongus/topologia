import numpy as np
import matplotlib.pyplot as plt


data = "estadistica.txt"
TotalSimu = 500
deep = 29

GenMean = []
StandDev = []

for col in range(deep):
    print("col: ", col)
    for i in range(TotalSimu):
        try:
            generations = np.loadtxt(data, usecols=col, max_rows=i)
        except:
            generations = np.loadtxt(data, usecols=col, max_rows=i - 1)
            print("Except:")
            break
    print(generations)
    GenMean.append(generations.sum())
    StandDev.append(generations.std(ddof=-(50 - generations.size)))
    print()

Statistic = np.array(GenMean) / TotalSimu
print(StandDev)

x = np.array([i for i in range(Statistic.size)])
plt.xlim(0, deep)
plt.xticks([2 * i for i in range(deep // 2 + 1)])
plt.ylim(0, 15)
plt.plot(x, Statistic, marker="o")
plt.grid()
plt.show()
