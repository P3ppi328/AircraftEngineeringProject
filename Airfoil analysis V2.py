# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 13:34:57 2023

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
import time
starttime = time.time()

        
read_polar_file = read_polar("NACA23118_val_V2")
foil = fileio.import_airfoil_data("NACA23118.txt")


# plt.figure()
# plt.plot(foil[0][0],foil[0][1], foil[1][0],foil[1][1])
# plt.ylim(-0.25,0.25)


foilx_lower = foil[1][0]
foily_lower = foil[1][1]

foilx_upper = foil[0][0]
foily_upper = foil[0][1]


# plt.figure()
# plt.plot(foil[0][0],foil[0][1], foil[1][0],foil[1][1])
# plt.ylim(-0.25,0.25)
# foil = Airfoil.NACA4('2412')

# foilx_lower = foil._x_lower
# foily_lower = foil._y_lower

# foilx_upper = foil._x_upper
# foily_upper = foil._y_upper



#constants
g = 9.81
alphacruise= 4
alphatake = 17
MASS = 74.67 #[kg]
rho_alu = 2710 #[kg/m^-3]  
LoD = 30 #L/D which can be specified
h = 3500 #[m] height of the aircraft at cruis (in most extreme conditions)
Cd= read_polar_file.get("cd")[463] #during cruise Cd must be lowest
#Cd= 0.007
#Cl = 0.9
Cl = read_polar_file.get("cl")[np.where(read_polar_file.get("cd") == Cd)[0][0]] # from XFoil for NACA2412 at alfa=2.5
Cl_max = 1.83
U = 150/3.6 #[m*s^-1] velocity of the plane at cruise
e = 0.8 # efficiency number wing 0.7 for rectengular 
Cr = 0.85
Tapered_ratio =0.3
Ct = Cr*Tapered_ratio




Chord_av = (Ct+Cr)/2

#properties
  
  #Air
AirProperties = Atmosphere(h)

#Mach & Reynolds number
Re = (AirProperties.density[0]*U*Chord_av)/AirProperties.dynamic_viscosity[0]
Ma = U/AirProperties.speed_of_sound[0]


#load cases at cruise
LIFT = (MASS*g)/2 #[N] during cruise



#calculate area (S) needed to generate W in L
WingA = LIFT/((AirProperties.density[0]/2)*(U**2)*Cl)

b = WingA/Chord_av


#Drag with wing area
D = Cd*((AirProperties.density[0]/2)*(U**2)*WingA)

AR = b**2/WingA

alphaicruise = Cl/(np.pi*AR)
alphaeffeccruise = alphacruise - alphaicruise

alphaitake = Cl_max/(np.pi*AR)
alphaeffectake = alphatake - alphaitake

D_ind = (Cl**2)/(np.pi *AR * e)*((AirProperties.density[0]/2)*(U**2)*WingA)

D_res2 = D + D_ind
#calculate possible AR to achieve set L/D taking Di as D_res - D
#AR = (Cl**2)/(np.pi *(D_res-D) * e)*((AirProperties.density[0]/2)*(U**2)*WingA)
#b = np.sqrt(AR*WingA)

# %% Determine the quarter-chord point and 2/3rd on the airfoil


        
# for i,_ in enumerate(foilx_lower):
#     if (abs(_-2*((LE+TE)/3)) <= 0.1):
#         index_rear_spar = i
   


# %% Chord function

b_span = np.linspace(0, b,len(foilx_lower)*2)
Chord_func = ((Ct-Cr)/b)*b_span+Cr





# # Make data
# foilx_lower = foil[1][0]
# foily_lower = foil[1][1]

# foilx_upper = foil[0][0]
# foily_upper = foil[0][1]


# plt.figure()
# plt.plot(foilx_lower,foily_lower)
# scaled = foilx_lower*Chord_func
# scaled2 = foily_lower*Chord_func
# tD_shapes_lower = np.array((scaled,scaled2))
# X,Y = np.meshgrid(foilx_lower,foily_lower)
# X,Z = np.meshgrid(foilx_lower,Chord_func)
# # Plot the surface
# ax.plot_surface(X,Y,Z)
# X,Y = np.meshgrid(foilx_upper,foily_upper)
# Z = Chord_func.reshape(X.shape)
# ax.plot_surface(X,Y,Z)
# ax.axes.set_xlim3d(left=-1, right=1) 
# ax.axes.set_ylim3d(bottom=-1, top=1) 
ax = plt.figure().add_subplot(projection='3d')
ax.view_init(elev=0, azim=-90, roll=-90)
ax.set_zlim3d(-4, 4)                     
ax.set_ylim3d(-2, 2)                     
ax.set_xlim3d(-2, 2) 
for i,_ in enumerate(Chord_func):
    



    # Plot a sin curve using the x and y axes.
    x = np.append(foilx_lower, foilx_upper)*_
    y = np.append(foily_lower, foily_upper)*_
    
    ax.plot(y,x , zs=b_span[i], zdir='z', c = "red",label='curve in (x, y)')
    
    

# Set an equal aspect ratio
#ax.set_aspect('equal')
Mass_front_spar = 0
Mass_skin = 0
Area_tot = np.array([])
thickness_airfoil_final = []
for i,chord in enumerate(Chord_func):
    
    
    #print(chord)

    # Plot a sin curve using the x and y axes.
    x = np.append(foilx_lower, foilx_upper)*chord
    y = np.append(foily_lower, foily_upper)*chord
    
    ax.plot(y, x, zs=-0.2-(b_span[i]), zdir='z', c = "red",label='curve in (x, y)')
    #Airfoil
    LE = x[0]
    TE = x[int(len(x)/2)]
    #quarter-chord point
    index_front_spar = np.array([])
    for index_chord4,_ in enumerate(x):
        if (abs(_-((LE+TE)/4)) <= 0.1):
            index_front_spar = np.append(index_front_spar,int(index_chord4))
        else:
            continue
            

    y_4cp_upper = y[int(index_front_spar[1])]
    y_4cp_lower = y[int(index_front_spar[0])]
    #thickness_airfoil = abs(y_4cp_upper)+abs(y_4cp_lower)
    thickness_airfoil = 0.18*chord
    print(thickness_airfoil)
    thickness_airfoil_final.append(thickness_airfoil)
    
    #print(thickness_airfoil)


    #2/3rds on the chord line
    # y_rear_spar_upper = y[index_rear_spar]
    # print(y_rear_spar_upper)
    # y_rear_spar_lower = foily_lower[index_rear_spar]
    # %% Material dimensions
    #spar
    thickness_spar = 10*10**-3
    width_front_spar = thickness_airfoil*chord
    
    if i < len(b_span)-1:
        dlength_spar = b_span[i+1]-b_span[i]
        dMass_front_spar = dlength_spar*width_front_spar*thickness_spar*rho_alu
        
        Mass_front_spar += dMass_front_spar
        #print(Mass_front_spar)
        #print(dlength_spar)
    else:
        continue
    #width spar
            
    
    
          #front spar
    #print(width_front_spar)
    # width_rear_spar = (y_rear_spar_upper-y_rear_spar_lower)*Chord #rear spar

    #ribs
    N_ribs = 10 #amount of ribs
    thickness_rib = 5*10**-3
    #Area of the airfoil
    Area_airfoil = 0
    for i,value in enumerate(foily_upper):
        if i == len(foily_upper)-1:
            break
        else:
            A = ((abs(foily_upper[i])+abs(foily_lower[i])))*((foilx_upper[i]-foilx_upper[i+1]))*chord**2
            #print((abs(y[i])+abs(y[-i])))
            Area_airfoil += A
    
    Area_tot = np.append(Area_tot,Area_airfoil)
    
    
    
    
        
    #print(Area_airfoil)
    #skin
    
    thickness_skin = 1*10**-3
    #Circumference Airfoil
    Cir = 0
    for i,_ in enumerate(foily_upper):
        if i == len(foily_upper)-1:
            break
        else:
            s = np.sqrt(((foily_upper[i+1]-foily_upper[i])*chord)**2+((foilx_upper[i+1]-foilx_upper[i])*chord)**2)+np.sqrt(((foily_lower[i+1]-foily_lower[i])*chord)**2+((foilx_upper[i+1]-foilx_upper[i])*chord)**2)
            Cir += s
    if i < len(b_span)-1:
        dlength_skin = b_span[i+1]-b_span[i]
        dMass_skin = Cir*dlength_skin*thickness_skin*rho_alu
        
        Mass_skin += dMass_skin
        #print(Mass_front_spar)
        #print(dlength_spar)
    else:
        continue
    
    
plt.show()



        
# %% Mass spars, ribs and skin
#spars

#Mass_rear_spar = length_spar*width_rear_spar*thickness_spar*rho_alu

#ribs
Mass_ribs = np.mean(Area_tot)*thickness_rib*N_ribs*rho_alu

#skin


#total mass
Total_wing_mass = Mass_front_spar+Mass_ribs+Mass_skin



MASS_0 = 74.67
diff_wing_mass = 0
Mass_front_spar_new = 0
b_span_list = np.array([])
while True:
    
    if MASS+diff_wing_mass==MASS_0:
        MASS+=Total_wing_mass
        
    else:
        
        
      
    
        #load cases at cruise
        LIFT = (MASS*g)/2 #[N] during cruise
        #print(LIFT/g)
        
    
        #calculate area (S) needed to generate W in L
        WingA_new = LIFT/((AirProperties.density[0]/2)*(U**2)*Cl)
        b_new = WingA_new/Chord_av
        
        
        AR_new = b_new**2/WingA_new
        #Drag with wing area
        D = Cd*((AirProperties.density[0]/2)*(U**2)*WingA_new)
        
        D_ind = (Cl**2)/(np.pi *AR_new * e)*((AirProperties.density[0]/2)*(U**2)*WingA_new)
        D_res = D + D_ind
        Cd_new_new=D_res/((AirProperties.density[0]/2)*(U**2)*WingA_new)
        #calculate possible AR to achieve set L/D taking Di as D_res - D
        
        
        
        diff_b = b_new - b
        
        b_span_new = np.linspace(0, b_new,len(foilx_lower)*2)
        
        Chord_func_new = ((Ct-Cr)/b_new)*b_span_new+Cr
        Mass_front_spar_new = 0
        Mass_skin_new = 0
        Area_tot_new = np.array([])
        for i,chord in enumerate(Chord_func_new):
            
            
            #print(chord)

            # Plot a sin curve using the x and y axes.
            x = np.append(foilx_lower, foilx_upper)*chord
            y = np.append(foily_lower, foily_upper)*chord
            
            ax.plot(y, x, zs=-0.2-(b_span[i]/2), zdir='z', c = "red",label='curve in (x, y)')
            #Airfoil
            LE = x[0]
            TE = x[int(len(x)/2)]
            #quarter-chord point
            index_front_spar = np.array([])
            for index_chord4,_ in enumerate(x):
                if (abs(_-((LE+TE)/4)) <= 0.1):
                    index_front_spar = np.append(index_front_spar,int(index_chord4))
                else:
                    continue
                    

            y_4cp_upper = y[int(index_front_spar[1])]
            y_4cp_lower = y[int(index_front_spar[0])]
            thickness_airfoil = abs(y_4cp_upper)+abs(y_4cp_lower)
       
            
            #print(thickness_airfoil)


            #2/3rds on the chord line
            # y_rear_spar_upper = y[index_rear_spar]
            # print(y_rear_spar_upper)
            # y_rear_spar_lower = foily_lower[index_rear_spar]
            # %% Material dimensions
            #spar
            thickness_spar = 10*10**-3
            width_front_spar = thickness_airfoil*chord
            
            if i < len(b_span_new)-1:
                dlength_spar_new = b_span_new[i+1]-b_span_new[i]
                dMass_front_spar_new = dlength_spar_new*width_front_spar*thickness_spar*rho_alu
                
                Mass_front_spar_new += dMass_front_spar_new
                #print(Mass_front_spar)
                #print(dlength_spar)
            else:
                continue
            #width spar
                    
            
            
                  #front spar
            #print(width_front_spar)
            # width_rear_spar = (y_rear_spar_upper-y_rear_spar_lower)*Chord #rear spar

            
                
            #print(Area_airfoil)
            #skin
            
            thickness_skin = 1*10**-3
            #Circumference Airfoil
            Cir = 0
            for i,_ in enumerate(foily_upper):
                if i == len(foily_upper)-1:
                    break
                else:
                    s = np.sqrt(((foily_upper[i+1]-foily_upper[i])*chord)**2+((foilx_upper[i+1]-foilx_upper[i])*chord)**2)+np.sqrt(((foily_lower[i+1]-foily_lower[i])*chord)**2+((foilx_upper[i+1]-foilx_upper[i])*chord)**2)
                    Cir += s
            if i < len(b_span)-1:
                dlength_skin_new = b_span_new[i+1]-b_span_new[i]
                dMass_skin = Cir*dlength_skin_new*thickness_skin*rho_alu
                
                Mass_skin_new += dMass_skin
                #print(Mass_front_spar)
                #print(dlength_spar)
            else:
                continue
            
            
       
        diff_mass_front_spar = Mass_front_spar_new-Mass_front_spar
        Mass_front_spar += diff_mass_front_spar
        
        diff_mass_skin = Mass_skin_new-Mass_skin
        Mass_skin += diff_mass_skin
        
        
        #
        diff_wing_mass = diff_mass_front_spar+diff_mass_skin
        #print(diff_wing_mass)
        
        b += diff_b
        Total_wing_mass += diff_wing_mass
        MASS += diff_wing_mass
        
        if diff_wing_mass<(0.0000001*Total_wing_mass):
            print(f"Current wing span (b): {b:.2f} m")
            print(f"Current total mass: {MASS:.2f} kg")
            print(f"Wing mass: {MASS-MASS_0:.2f} kg\n")
            print(f'Drag: {np.round(D_res,2)} N')
            print(f'LoD: {np.round((MASS*9.81)/(2*D_res),2)}')
            print(f'runtime = {np.round(time.time() - starttime,3)} s')
            break

   
    