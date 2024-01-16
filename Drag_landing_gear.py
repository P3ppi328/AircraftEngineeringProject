# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 10:45:12 2024

@author: pepij
"""



from airfoils import fileio
from airfoils import Airfoil
from ambiance import Atmosphere
import numpy as np
from NACA5_generator import naca5
import matplotlib.pyplot as plt
import materials
from xfoil_reading import read_polar
from scipy.integrate import simpson
#constants
g = 9.81
alphacruise= 4
alphatake = 17
MASS = 80 #[kg]
h = 3500 #[m] height of the aircraft at cruis (in most extreme conditions)
U = 150/3.6 #[m*s^-1] velocity of the plane at cruise
AirProperties = Atmosphere(h)
Cr = 0.7
Tapered_ratio =0.3
Ct = Cr*Tapered_ratio
Cd = 0.00773

d_rod = 0.075
l=0.175
chord_rod = d_rod/0.2


foil = fileio.import_airfoil_data("NACA0020.txt")
foilx_lower = foil[1][0]
foily_lower = foil[1][1]

foilx_upper = foil[0][0]
foily_upper = foil[0][1]


#Circumference Airfoil
Cir = 0
for i,_ in enumerate(foily_upper):
    if i == len(foily_upper)-1:
        break
    else:
        s = np.sqrt(((foily_upper[i+1]-foily_upper[i])*chord_rod)**2+((foilx_upper[i+1]-foilx_upper[i])*chord_rod)**2)+np.sqrt(((foily_lower[i+1]-foily_lower[i])*chord_rod)**2+((foilx_upper[i+1]-foilx_upper[i])*chord_rod)**2)
        Cir += s

WingA=Cir*l

D_rod = Cd*((AirProperties.density[0]/2)*(U**2)*WingA)




