from airfoils import fileio
from airfoils import Airfoil
from ambiance import Atmosphere
import numpy as np
from NACA5_generator import naca5
import matplotlib.pyplot as plt
import materials
from xfoil_reading import read_polar
from scipy.integrate import simpson
import time

h = 250
Ucruise = 150/3.6
rho = Atmosphere(h).density[0]
g = Atmosphere(h).grav_accel[0]
U = np.linspace(0,Ucruise,100000)
Clmax = 2.23
WingA = 1.2461514265700173*2
MASS = 90
Ud = 175/3.6
W = MASS*g
Vstall = np.sqrt((2*W)/(rho*WingA*Clmax))
Lmaxdive = (rho*Ud**2*Clmax*WingA)/2
Lmaxcruise = (rho*Ucruise**2*Clmax*WingA)/2
limit_load_factor = 2
n2 = 0.75*limit_load_factor
R = Ucruise**2/(g*np.sqrt(limit_load_factor**2-1))



L = np.array([])
for i,_ in enumerate(U):
    if _ <= Ucruise and ((rho*_**2*Clmax*WingA)/2)/W < limit_load_factor:
        L = np.append(L,(rho*_**2*Clmax*WingA)/2)
    else:
        L = np.append(L,L[i-1])


Ln = np.array([])
for i,_ in enumerate(U):
    if _ < Vstall:
        Ln = np.append(Ln,-(rho*_**2*Clmax*WingA)/2)
    else:
        Ln = np.append(Ln, Ln[i-1])
        
a = (n2-L[-1]/W)/(Ud-Ucruise)
U2 = np.linspace(Ucruise,Ud,1000)
b = 1/(Ud-Ucruise)
Lf = (b*(U2-Ucruise)-1)
Lff = (a*(U2-Ucruise)+limit_load_factor)
plt.figure()
plt.plot(U,L/W, color="red")
plt.plot(U,Ln/W, color="red")
plt.vlines(Vstall,-1,1, label="Stall speed", color="blue")
# plt.hlines(limit_load_factor, 200,250)
plt.hlines(1,Vstall,Ud, linestyle="dashed", label="Level flight", color="purple")
plt.vlines(Ud,0,n2, color="red")
plt.vlines(Ucruise, -1,limit_load_factor, linestyle="dashed",label = "Cruise speed", color="green")
plt.plot(U2,Lf, color="red")
plt.plot(U2, Lff, color="red")
plt.xlabel("Speed [$m*s^{-1}$]")
plt.ylabel("load factor (n) [arb. units]")
plt.legend()

plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
