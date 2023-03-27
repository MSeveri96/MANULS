#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
import matplotlib




font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   :  14}

matplotlib.rc('font', **font)

######################################################################
##################### USER DEFINED SECTION  ##########################
######################################################################

# This is the only section that needs to be modified, the rest is up to the code

# Loading of the files containing the data, please see the documentation
data=np.loadtxt('/home/marco/cumulene/scan_polar_long_grid.txt')
polar=np.loadtxt('/home/marco/cumulene/scan_polar_long_polarizability.txt')


# "pointsx" and "pointsy" = number of points scanned for each variable
pointsx=int(37)

#the optimal electric field
#e=np.array([field_x,  field_y,  field_z])
e=np.array([-0.00395539,  0.00419754,  0.00640764])




######################################################################
##################### END OF USER DEFINED SECTION  ###################
######################################################################




########################################################################
########################################################################
########################################################################



energy=data[:,1]
dipx=data[:,2]
dipy=data[:,3]
dipz=data[:,4]

energy_column=data[:,1]


a00=np.array([])
a01=np.array([])
a02=np.array([])
a10=np.array([])
a11=np.array([])
a12=np.array([])
a20=np.array([])
a21=np.array([])
a22=np.array([])
for i in range(int(np.shape(polar)[0]/3)):
    a00=np.append(a00,polar[0+i*3,0])
    a01=np.append(a01,polar[0+i*3,1])
    a02=np.append(a02,polar[0+i*3,2])
    a10=np.append(a10,polar[1+i*3,0])
    a11=np.append(a11,polar[1+i*3,1])
    a12=np.append(a12,polar[1+i*3,2])
    a20=np.append(a20,polar[2+i*3,0])
    a21=np.append(a21,polar[2+i*3,1])
    a22=np.append(a22,polar[2+i*3,2])
    
    


fig=plt.figure()
x=np.linspace(data[0,0],data[-1,0],num=pointsx)
plt.plot(x,(energy-np.amin(energy))*627.503,marker='o')
plt.title('Unperturbed PES')
plt.xlabel(r'Variable Scanned')
plt.ylabel(r'Relative Energy (kcal/mol) ')
plt.show()
         
            
        
            
        

#####POLARIZABILITY
        
energy_perturbed=np.zeros(pointsx)





for g in range(int(pointsx)):
        energy_perturbed[g]=energy[g]-(dipx[g]*e[0]+dipy[g]*e[1]+dipz[g]*e[2])-(1/2)*e[0]*(e[0]*a00[g]+e[1]*a01[g]+e[2]*a02[g])-(1/2)*e[1]*(e[0]*a10[g]+e[1]*a11[g]+e[2]*a12[g])-(1/2)*e[2]*(e[0]*a20[g]+e[1]*a21[g]+e[2]*a22[g])
fig=plt.figure()
x=np.linspace(data[0,0],data[-1,0],num=pointsx)
plt.plot(x,(energy_perturbed-np.amin(energy))*627.503,marker='o',color='tab:orange')
plt.title('Perturbed PES')
plt.xlabel(r'Variable Scanned')
plt.ylabel(r'Relative Energy (kcal/mol)')
plt.show()
            
            
