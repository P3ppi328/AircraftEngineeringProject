# code propulsion system comparison

# Author: Minze Hamstra S2569620

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ambiance import Atmosphere

# propeller 1: Tmax = 23.27   https://shop.mejzlik.eu/folding-propeller-32x10-6-ccw-2b-mc-puller/
# propeller 2: Tmax = 9.35    https://shop.mejzlik.eu/folding-propeller-22x7-4-ccw-2b-mc-puller/
# propeller 3: Tmax = 3.51    https://shop.mejzlik.eu/folding-propeller-18x6-ccw-2b-mc-puller/      3 propellors needed
# propeller 4: Tmax = 5.72    https://www.vertiq.co/vert8108-220kv-mej22x7-4-24v                    2 propellors needed

T =  35.2 # [kgf] propeller 1
n = 1 # number of propellers

max_alpha = True

h_0 = 0

ambiance = Atmosphere(h_0)
rho = ambiance.density[0]
vis = ambiance.dynamic_viscosity[0]
g = ambiance.grav_accel[0]

# plane characteristics
M_total = 102                # [kg] total mass
M_payload = 25              # [kg] payload mass

h_wing = h_0 +1.2             # [m] vertical distance between the plane and the ground at take off
span = 2.26*2               # [m] wingspan
S_wing = 2*1.25             # [m^2] Wing surface 2D element
AR = 8.16                   # [] aspect ratio of the main wing 12/12/2023
e = 0.8                     # [] Oswalt efficiency number (0.75-0.8 for fuselage, 0.9-0.95 for main wing)

Cl = 0.4694                 # [] lift coefficient
Cl_max = 2.23               # [] max lift coefficient, this is of the 2d wing element. Cl max of the whole plane should be a bit smaller

V_cruise = 150/3.6          # [m/s] Cruise velocity of the plane
V_stall = ((2*M_total*g) / (rho*S_wing*Cl_max))**0.5             # [m/s] stall velocity of the plane, assumption and is not confirmed by aero
V_stall_empty = ((2*(M_total-M_payload)*g) / (rho*S_wing*Cl_max))**0.5             # [m/s] stall velocity of the plane, assumption and is not confirmed by aero
V_takeoff = 1.2*V_stall     # [m/s] take off velocity
V_landing = 1.3*V_stall

print('Vstall = ',V_stall)
print('V takeoff = ',V_takeoff)

l_fus = 3                   # [m] fuselage length
S_fugelage = 0.25*l_fus      # [m^2] Fugelage surface, approximated by a cylinder 06/12/2023

Fd_tail = 13.4                                                  # [N] tail drag given by Martin
# Drag coefficient tail
AR_HT = 4                   # horizontal AR
AR_VT=2                     # vertical AR
S_vert_tail = 0.43
S_hori_tail = 2*0.51

it = 2*(np.pi/180)          # inclination 
at =6.3                     # slope lift curve
Cl_tail = it*at             # lift coefficient tail
Cd0_tail = 0.01             # zero lift drag coefficient tail airfoil
Cdi_tail = (Cl_tail**2)/(np.pi*AR_HT*e); # induced drag coefficient horizontal tail section
Cd_tail_horizontal = Cdi_tail+Cd0_tail # drag coefficient horizontal tail section


Cd_landing_gear = 0.00773
S_landing_gear_rear = 0.011
S_landing_gear_front = 2 * 0.3855
S_wheels_front = 0.14
S_landing = S_landing_gear_front+ S_landing_gear_rear+ S_wheels_front
mu_tire = 0.05

h_altitude = 3500

Fthrust = T * n * g


def drag(h,v,Cl,Cl_max,l_fus,S_fugelage,AR,S_wing,e,max_alpha,Cd0_tail,Cd_tail_horizontal,S_vert_tail,S_hori_tail,Cd_landing_gear,S_landing,mu_tire,M_total):
    ambiance = Atmosphere(h)
    rho = ambiance.density[0]
    vis = ambiance.dynamic_viscosity[0]
    g = ambiance.grav_accel[0]
    h_wing = h + 1.2

    if max_alpha:
        Cl_value = Cl_max
    else:
        Cl_value = Cl

    GE_cof = (16*h_wing/span)**2 / (1 + (16*h_wing/span)**2)        # [] ground effect coefficient
    Cd_2d = 0.00857                                                 # [] drag coefficient 2D wing element,(zero lift) aproximation of aero 09/01/2024
    Cd_3d = GE_cof * (Cl_value**2)/(np.pi*AR*e)                     # [] drag coefficient lift induced drag lift off
    Cd_wing = Cd_2d + Cd_3d

    Re = (rho * V_cruise*l_fus)/vis
    Cd_fuselage = 0.074/ (Re**(1/5))                                # fuselage drag coefficient

    Fd_wing = Cd_wing* 0.5 * rho * v**2 * S_wing 
    Fd_fuselage = Cd_fuselage * 0.5 * rho * v**2 * S_fugelage

    # print()
    # print('Cd_2d main wing = ',Cd_2d)
    # print('Cd 3d wing = ' , Cd_3d)
    # print('Cd wing total = ', Cd_wing)
    # print('Cd fuselage = ',Cd_fuselage)
    # print('Cd 2d tail = ',Cd0_tail)
    # print('Cd 3d tail = ', Cdi_tail)
    # print('Cd horizontal tail = ',Cd_tail_horizontal)
    # print('Cd landing gear = ',Cd_landing_gear)
    # print('mu rolling = ', mu_tire)
    # print()

    Fd_tail = 13.4
    Fd_tail_vert = Cd0_tail * 0.5 * rho * v **2 * S_vert_tail
    Fd_tail_hori = Cd_tail_horizontal * 0.5 * rho * v **2 * S_hori_tail
    Fd_tail = Fd_tail_vert + Fd_tail_hori

   
    Fd_landing_gear = Cd_landing_gear * 0.5 * rho * v**2 * S_landing                                          # [N] landing gear drag force given by JP

    if h ==0 : 
        Fd_rolling =  0.5*M_total*g*mu_tire 
    else:
        Fd_rolling =0

    Fd = Fd_wing + Fd_fuselage + Fd_tail + Fd_landing_gear + Fd_rolling
    
    print('Fd wing = ', Fd_wing)
    print('Fd fuselage = ',Fd_fuselage)
    print('Fd tail = ',Fd_tail)
    print('Fd landing gear = ',Fd_landing_gear)
    print('Fd rolling = ',Fd_rolling)
    print('Fd = ',Fd)
    print()
    print()
    return Fd

def s_takeoff(V_takeoff,M_total,F_propeller_low,Fd):
    s = (V_takeoff**2 * M_total)/(2*(F_propeller_low - Fd))
    return s

def lift(h,v,Cl,S_wing):
    ambiance = Atmosphere(h)
    rho = ambiance.density[0]
    L = Cl * 0.5 * rho * v **2 * S_wing
    return L





print('takeoff')
Fd_takeoff = drag(h_0,V_takeoff,Cl,Cl_max,l_fus,S_fugelage,AR,S_wing,e,True,Cd0_tail,Cd_tail_horizontal,S_vert_tail,S_hori_tail,Cd_landing_gear,S_landing,mu_tire,M_total)
print('cruise')
Fd_altitude = drag(250,(200/3.6),Cl,Cl_max,l_fus,S_fugelage,AR,S_wing,e,False,Cd0_tail,Cd_tail_horizontal,S_vert_tail,S_hori_tail,Cd_landing_gear,S_landing,mu_tire,M_total)
print()
print('drag takeoff = ',Fd_takeoff)
print('drag cruise = ',Fd_altitude)


s = s_takeoff(V_takeoff,M_total,Fthrust,Fd_takeoff)
print()
print('takeoff distance = ',s)

L = lift(h_0,V_takeoff,Cl_max,S_wing)
W = g * M_total
print()
print('lift = ',L)
print('Weight = ',W)

Fd_landing = drag(0,V_landing,Cl,Cl_max,l_fus,S_fugelage,AR,S_wing,e,True,Cd0_tail,Cd_tail_horizontal,S_vert_tail,S_hori_tail,Cd_landing_gear,S_landing,mu_tire,M_total)

s_landing = (V_landing**2 *M_total) / (2*Fd_landing) 
print()
print('Landing distance = ',s_landing)


s_given = 150
F = 0.5* ( V_takeoff**2 *M_total)/s_given + Fd_takeoff
print('thrust required takeoff = ',F)
