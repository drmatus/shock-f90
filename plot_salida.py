#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import glob

files = glob.glob("out*.dat")
files.sort()

n = 30
i = 0

#anhco y alto del plot  (en pulgadas)
ancho = 12.0/2.54
alto = 17/2.54

for i,nombre in enumerate(files):

    print("Reading ", nombre)

    #Leer el tiempo del snapshot
    data = open(nombre,"r")
    # for j in range(2):
    #     tmp = data.readline()
    header = data.readline()
    t = float(header.split()[1])
    data.close()


    #Leer los datos del snapshot
    x, rho, v, P = np.genfromtxt(nombre, unpack=True, skip_header=3)

    nombre = "graficos_{:06d}.png".format(i)
    print(nombre)

    fig, (ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, dpi=300,
            figsize=(ancho,alto))

    titulo = "t = {:f} yr".format(t/3.1e7)
    fig.suptitle(titulo)
    ax1.plot(x,v)
    ax1.set_ylabel("Velocidad [cm s-1]")

    ax2.plot(x,P)
    ax2.set_ylabel("Presion [erg cm${^-3}$]")
    ax2.set_yscale('log')

    ax3.plot(x,rho)
    ax3.set_ylabel("Densidad [gr cm${^-3}$]")
    ax3.set_yscale('log')

    ax3.set_xlabel("Postion [cm]")

    fig.subplots_adjust(left=0.20, right=0.99,top=0.95,hspace=0.0)

    plt.savefig(nombre)
    plt.close(fig)

    # if(i > 50):
    #     break

