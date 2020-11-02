import matplotlib.pyplot as plt
import numpy as np
import sys

file = open("int/pyinp.txt", "r")

action = file.readline()

dd = True
pp = True

if(action == "n\n"):
    dd = False
    pp = False
if(action == "d\n"):
    pp = False
if(action == "p\n"):
    dd = False

if(dd or pp):    
    moCount = int(file.readline())
    xStart = float(file.readline())
    yStart = float(file.readline())
    xNum = int(file.readline())
    yNum = int(file.readline())
    step = float(file.readline())
    name = file.readline()
    name = name[:len(name)-1]

    xEnd = xStart + step * (xNum-1)
    yEnd = yStart + step * (yNum-1)
    Num = 0

    if(xNum > yNum):
        yEnd = yStart + step * (xNum-1)
        Num = xNum
    else:
        xEnd = xStart + step * (yNum-1)
        Num = yNum
    file.close()
    if(pp):
        for a in range(0, moCount):
            y, x = np.meshgrid(np.linspace(xStart, xEnd , Num), np.linspace(yStart, yEnd, Num))
            z = x*y
            f = open("int/plain_"+ str(a) +".txt", "r")
            for i in range(0, Num):
                for j in range(0, Num):
                    if(i >= xNum or j >= yNum):
                        z[i,j] = 0
                    else:
                        z[i,j] = f.readline()
            f.close()
            z = z[:-1, :-1]
            z_min, z_max = -np.abs(z).max(), np.abs(z).max()

            fig, ax = plt.subplots()

            c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
            ax.set_title(name+"_"+str(a))

            ax.axis([x.min(), x.max(), y.min(), y.max()])
            ax.set_ylabel('Distance [A]')
            ax.set_xlabel('Distance [A]')
            fig.colorbar(c, ax=ax)
            plt.savefig("pic/"+name+"_"+str(a)+".png")
            plt.clf()

    if(dd):
        y, x = np.meshgrid(np.linspace(xStart, xEnd , Num), np.linspace(yStart, yEnd, Num))
        z = x*y
        f = open("int/dplain.txt", "r")
        for i in range(0, Num):
            for j in range(0, Num):
                if(i >= xNum or j >= yNum):
                    z[i,j] = 0
                else:
                    z[i,j] = f.readline()
        f.close()
        z = z[:-1, :-1]
        z_min, z_max = -np.abs(z).max(), np.abs(z).max()

        fig, ax = plt.subplots()

        c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
        ax.set_title(name+"_density")

        ax.axis([x.min(), x.max(), y.min(), y.max()])
        ax.set_ylabel('Distance [A]')
        ax.set_xlabel('Distance [A]')
        fig.colorbar(c, ax=ax)
        plt.savefig("pic/"+name+"_density"+".png")
        plt.clf()
else:
    file.close()


