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

data=np.loadtxt('/path/to/file_grid.txt')
polar=np.loadtxt('/path/to/file_polarizability.txt')



# "pointsx" and "pointsy" = number of points scanned for each variable
pointsx=int(23)
pointsy=int(23)

#the optimal electric field
#e=np.array([field_x,  field_y,  field_z])
e=np.array([-0.00395539,  0.00419754,  0.00640764])




######################################################################
##################### END OF USER DEFINED SECTION  ###################
######################################################################




########################################################################
########################################################################
########################################################################



energy=np.reshape(data[:,2],(pointsy,pointsx),order='F')

dipx=np.reshape(data[:,3],(pointsy,pointsx),order='F')
dipy=np.reshape(data[:,4],(pointsy,pointsx),order='F')
dipz=np.reshape(data[:,5],(pointsy,pointsx),order='F')

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
    

a00=np.reshape(a00,(pointsy,pointsx),order='F')
a01=np.reshape(a01,(pointsy,pointsx),order='F')
a02=np.reshape(a02,(pointsy,pointsx),order='F')
a10=np.reshape(a10,(pointsy,pointsx),order='F')
a11=np.reshape(a11,(pointsy,pointsx),order='F')
a12=np.reshape(a12,(pointsy,pointsx),order='F')
a20=np.reshape(a20,(pointsy,pointsx),order='F')
a21=np.reshape(a21,(pointsy,pointsx),order='F')
a22=np.reshape(a22,(pointsy,pointsx),order='F')


fig=plt.figure()
x=np.linspace(data[0,0],data[-1,0],num=pointsx)
y=np.linspace(data[0,1],data[-1,1],num=pointsy)
levels = 30
cmap= plt.cm.get_cmap("viridis", levels+1)
plt.contour(x,y,(energy-np.amin(energy))*627.503,levels,cmap=cmap)
plt.title('Unperturbed PES')
plt.xlabel(r'Variable 1 Scanned')
plt.ylabel(r'Variable 2 Scanned')
cbar=plt.colorbar()
cbar.set_label('Relative Energy (kcal/mol)')
plt.show()
         
            
        
            
        

#####POLARIZABILITY
        
energy_perturbed=np.zeros((pointsy,pointsx))




print(e)
for m in range(int(pointsy)):
    for n in range(int(pointsx)):
        energy_perturbed[m,n]=energy[m,n]-(dipx[m,n]*e[0]+dipy[m,n]*e[1]+dipz[m,n]*e[2])-(1/2)*e[0]*(e[0]*a00[m,n]+e[1]*a01[m,n]+e[2]*a02[m,n])-(1/2)*e[1]*(e[0]*a10[m,n]+e[1]*a11[m,n]+e[2]*a12[m,n])-(1/2)*e[2]*(e[0]*a20[m,n]+e[1]*a21[m,n]+e[2]*a22[m,n])


fig=plt.figure()
x=np.linspace(data[0,0],data[-1,0],num=pointsx)
y=np.linspace(data[0,1],data[-1,1],num=pointsy)
levels = 30
cmap= plt.cm.get_cmap("viridis", levels+1)
plt.contour(x,y,(energy_perturbed-np.amin(energy))*627.503,levels,cmap=cmap)
plt.title('Perturbed PES')
plt.xlabel(r'Variable 1 Scanned')
plt.ylabel(r'Variable 2 Scanned')
cbar=plt.colorbar()
cbar.set_label('Relative Energy (kcal/mol)')
plt.show()
