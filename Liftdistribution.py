import numpy as np
import matplotlib.pyplot as plt
import vortex 
from aeropy import aero_module
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

h=3500
U = 150/3.6

AirProperties = Atmosphere(h)
N=100
g = 9.81
tap4= aero_module.LLT_calculator(-1.1,0.0089,N=N, b=2.47, taper=0.4, chord_root=1, alpha_root=0.0, V=1.0)
tap3 = aero_module.LLT_calculator(-1.1,0.0089,N=N, b=2.47, taper=0.3, chord_root=1, alpha_root=0.0, V=1.0)
# tap0= aero_module.LLT_calculator(-1.1,0.0089,N=N, b=2.47, taper=1, chord_root=1, alpha_root=0.0, V=1.0)
dis4 = tap4.get("distribution")
dis3 = tap3.get("distribution")
# dis0 = tap0.get("distribution")
# np.put(dis,[-1],0.3)
Cl = tap4.get("cls")
b = 2.26
Cr = 0.85
Tapered_ratio =0.3

Ct = Cr*Tapered_ratio
b_span = np.linspace(0, b,N)
Chord_func = ((Ct-Cr)/b)*b_span+Cr
LIFT = 92*g
tot = sum(dis3)
dis44 = (dis4/tot)*LIFT
M = dis44*b_span
force = (dis3/tot)*LIFT
# LIFT= ((AirProperties.density[0]/2)*(U**2)*Cl*Chord_func)
plt.figure()
plt.plot(b_span, force , color="red")
plt.figure()
plt.plot(b_span,M)
total = np.transpose(np.vstack((b_span,force)))
#total = np.concatenate((b_span,force),axis=1)
# plt.plot(b_span,dis3, color="blue")
# plt.plot(b_span,dis0, color="purple")
# plt.ylim(0,2.5)
#from aero_module import LLT_calculator

#wing = aeropy.aero_module.LLT_calculator(alpha_L_0_root, c_D_xfoil, N=10, b=10.0, taper=1.0, chord_root=1, alpha_root=0.0, V=1.0)
#wing = aeropy.aero_module.LLT_calculator()
# def calculate_lift_distribution(span, taper_ratio, lift_total, n_panels=20):
#     # Create tapered wing
#     wing = vortex. .create_tapered_wing(span, taper_ratio, n_panels)

#     # Solve for the lift distribution
#     lift_distribution = lifting_line.solve(wing, lift_total)

#     return wing, lift_distribution

# def plot_lift_distribution(wing, lift_distribution):
#     # Plotting
#     plt.figure(figsize=(10, 6))
#     plt.plot(wing.y, lift_distribution, label='Lift Distribution')
#     plt.xlabel('Spanwise Distance (m)')
#     plt.ylabel('Lift per unit span (N/m)')
#     plt.title('Lift Distribution over a Tapered Wing')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# if __name__ == "__main__":
#     # Given parameters
#     span = 10  # Wing span in meters
#     taper_ratio = 0.4  # Taper ratio
#     lift_total = 1500  # Total lift in Newtons

#     # Calculate lift distribution
#     wing, lift_distribution = calculate_lift_distribution(span, taper_ratio, lift_total)

#     # Plot the lift distribution
#     plot_lift_distribution(wing, lift_distribution)
