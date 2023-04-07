# Grafica la epidemia en el tiempo.

import numpy as np
import matplotlib.pyplot as plt

data = ["./data/epid.txt"]

maxrows = None
N = 200

healthy = np.loadtxt(data[0], usecols=0, max_rows=maxrows)
infected = np.loadtxt(data[0], usecols=1, max_rows=maxrows)
refractary = np.loadtxt(data[0], usecols=2, max_rows=maxrows)
time = np.loadtxt(data[0], usecols=3, max_rows=maxrows)

plt.xlabel("Time")
plt.ylabel("Populations")
plt.ylim(0, N)

plt.axhline(y=0, color="black")
plt.axvline(x=0, color="black")

plt.plot(time, healthy, label="healthy")
plt.plot(time, infected, label="Infected")
plt.plot(time, refractary, label="Refractary")

plt.legend()
plt.grid()
plt.show()
