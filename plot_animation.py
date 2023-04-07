# Animación de la epidemia a través de agentes.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches

colores = [
    "blue",
    "red",
    "green",
]
archivo = "data/animacion.txt"

##############################################################################################
# Animacion
##############################################################################################


def trayectoria(np_steps, tpause=0.01):

    N = 200
    L = 70

    fig, ax = plt.subplots()

    for i in range(np_steps):

        if i % 100 == 0:
            print(i)

        x = np.loadtxt(archivo, usecols=0, skiprows=N * i, max_rows=N)
        y = np.loadtxt(archivo, usecols=1, skiprows=N * i, max_rows=N)

        estado = np.loadtxt(archivo, usecols=3, skiprows=N * i, max_rows=N, dtype=int)

        plt.cla()

        plt.title("Agents system")
        plt.xlabel("x coordinate")
        plt.ylabel("y coordinate")

        plt.axis("square")
        plt.grid()
        plt.xlim(-1, L + 1)
        plt.ylim(-1, L + 1)

        for j in range(N):
            circ = patches.Circle((x[j], y[j]), 1, alpha=0.7, fc=colores[estado[j]])
            ax.add_patch(circ)

        plt.savefig("video/pic%.4i.png" % (i), dpi=150)
        # plt.pause(tpause)


trayectoria(4000)


####################################################################
####################################################################
