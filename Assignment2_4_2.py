# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:50:01 2023

@author: pepij
"""

from airfoils import fileio
import numpy as np
import pandas as df
from xfoil_reading import read_polar
import matplotlib.pyplot as plt
# foil = fileio.import_airfoil_data("NACA632415.txt")

read_polar_file = read_polar("NACA23018_val_V3")
read_polar_file_V2 = read_polar("NACA23118_val_V2")
# foil_upper = np.flip(foil[0])
# foil_lower = foil[1]
# b = np.transpose(foil_upper)
# # c = np.array(b[0])

cl = read_polar_file.get("cl")
alpha = read_polar_file.get("a")

cl_V2 = read_polar_file_V2.get("cl")
alpha_V2 = read_polar_file_V2.get("a")

plt.figure()
plt.scatter(alpha,cl, label=r"$\alpha$ vs $C_{L} NACA23018$",s=0.5)
plt.scatter(alpha_V2,cl_V2, label=r"$\alpha$ vs $C_{L} NACA23118$",s=0.5, color='red')
plt.xlim(-25,25)
plt.ylim(-1.5,2.4)
plt.title(r"$\alpha$ vs $C_{L}$ for the NACA23118")
plt.ylabel(r"$C_{L}$")
plt.xlabel(r"$\alpha$ [arb. units]")
plt.scatter([16.9],[2.23], s=20,color="green", edgecolor='black',label=r'$C_{L,max} = 2.23$ with a fixed slat')
plt.vlines(alpha[np.where(cl==max(cl))[0][0]],max(cl), 2.23,color='green',linestyles='dashdot')
plt.vlines(alpha[np.where(cl==max(cl))[0][0]],-1.5,max(cl),color='purple',linestyles='dashed',label=f"cl max = {np.round(max(cl),2)} at AoA = {np.round(alpha[np.where(cl == max(cl))][0],2)} [arb. units]")
# plt.vlines(x[-1],0,max(v), colors="red", linestyles="dashed", label=r"$\nu_{max}$" f" = {np.round(max(v),3)} m")
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.

plt.minorticks_on()
plt.legend()

cd = read_polar_file.get("cd")

cd_slat = cd*1.5

plt.figure()
plt.scatter(alpha,cd, label=r"$\alpha$ vs $C_{D} NACA23018$",s=0.5)
plt.scatter(alpha,cd_slat, label=r"$\alpha$ vs $C_{D} NACA23118$",s=0.5, color='red')
plt.xlim(-25,25)
plt.ylim(0,0.27)
plt.title(r"$\alpha$ vs $C_{D}$ for the NACA23118")
plt.ylabel(r"$C_{D}$")
plt.xlabel(r"$\alpha$ [arb. units]")
#plt.scatter([16.9],[1.826], s=20,color="green", edgecolor='black',label=r'$C_{L,max} = 1.83$ with a fixed slat')
#plt.vlines(alpha[np.where(cl==max(cl))[0][0]],max(cl), 1.826,color='green',linestyles='dashdot')
#plt.vlines(alpha[np.where(cl==max(cl))[0][0]],-1.5,max(cl),color='purple',linestyles='dashed',label=f"cl max = {np.round(max(cl),2)} at AoA = {np.round(alpha[np.where(cl == max(cl))][0],2)} [arb. units]")
# plt.vlines(x[-1],0,max(v), colors="red", linestyles="dashed", label=r"$\nu_{max}$" f" = {np.round(max(v),3)} m")
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.

plt.minorticks_on()
plt.legend()


cm = read_polar_file.get("cm")


cm_V2 = read_polar_file_V2.get("cm")


plt.figure()
plt.scatter(alpha,cm, label=r"$\alpha$ vs $C_{M} NACA23018$",s=0.5)
plt.scatter(alpha_V2,cm_V2, label=r"$\alpha$ vs $C_{M} NACA23118$",s=0.5, color='red')
plt.title(r"$\alpha$ vs $C_{M}$ for the NACA23018 and the NACA23118")
plt.ylabel(r"$C_{M}$ [arb. units]")
plt.xlabel(r"$\alpha$ [arb. units]")
plt.xlim(-25,25)
plt.hlines(cm[39],-25,0,color='green',linestyles='dashdot',label=r'$C_{M}|_{AoA = 0}$' f' for the NACA23018 = {np.round(cm[39],4)}')
plt.hlines(cm_V2[386],-25,0,color='purple',linestyles='dashed',label=r'$C_{M}|_{AoA = 0}$' f'for the NACA23118  = {np.round(cm_V2[386],4)}')
# plt.vlines(x[-1],0,max(v), colors="red", linestyles="dashed", label=r"$\nu_{max}$" f" = {np.round(max(v),3)} m")
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)

# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
plt.legend()
