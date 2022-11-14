#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 08:48:58 2022

@author: Marco Severi
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy as sp
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   :  14}

matplotlib.rc('font', **font)


# Given the construction of the code it is very important to scan from the shorter
# distance to the larger distance (i.e scan from 1.0 Angstroms to 2.0 Angstroms and
# not vice versa)

data=np.loadtxt('scan_final_grid.txt')
# data=np.flip(data,axis=0)


# in "step_x" and "step_y" I store the size of the step in each direction

step_x=0.005
step_y=0.005

#  in "pointsx" and "pointsy" I store the number of points scanned for each variable
pointsx=int(9)
pointsy=int(5)

# Starting from the data written in columns I rearrange them in matrices, in 
# order to implement easily the derivatives. Also in this case the ordering
# "Fortran style" is very important. At the end the columns reprensent
#  the x-axis, while the rows represent the y-axis

energy=np.reshape(data[:,2],(pointsy,pointsx),order='F')
dipx=np.reshape(data[:,3],(pointsy,pointsx),order='F')
dipy=np.reshape(data[:,4],(pointsy,pointsx),order='F')
dipz=np.reshape(data[:,5],(pointsy,pointsx),order='F')




# In this case we have a small grid that does not contain the TS and OP1 (OP1 is
# the intermediate) points therefore I compute the energy of those points 
# separately using the geometries reported in J. AM. CHEM. SOC. 2006, 128, 2394-2409 
energy_ts=-614.652371082026
energy_op1=-614.685780937115

# I scale the energies with respect to OP1
energy=energy-energy_op1

energy_ts=energy_ts-energy_op1

# I store the energies also in column format, it will be handy later
energy_column=data[:,2]
energy_column=energy_column-energy_op1

energy_op1=0


# PLOTTTING THE PES IN 2D

fig=plt.figure()
x=np.linspace(data[0,0],data[-1,0],num=pointsx)
y=np.linspace(data[0,1],data[-1,1],num=pointsy)

levels = 50
cmap= plt.cm.get_cmap("viridis", levels+1)

plt.contour(x,y,energy*627.503, levels,cmap=cmap)

cbar=plt.colorbar()
cbar.set_label('Relative Energy (kcal/mol)')
plt.xlabel(r'P-C')
plt.ylabel(r'C-O')





# in "pertubed_gradient" I store the gradient of the pes modified by the field
# this should be zero for the OBBP
perturbed_gradient=np.zeros(2)

# in "fields" I store the electric fields calculated for each point of the grid
fields=np.zeros(3)

# in "gradient_pes" I store the gradient of the unperturbed pes
gradient_pes=np.zeros(2)

# in "H_e" I store the Hessian matrices modified by the electric field.
# It is (H-Di*ei*E)
H_e=np.zeros((1,2,2))

# in "indexes" I store the indexes used for looping through the points
# of the pes
indexes=np.zeros(2)

# in "norm_s" i store the values of the norm of the vector "s" (eq 22 ref. work)
norm_s=np.array([])

# in "s" i store the values of the vector "s" (eq 22 ref. work)
s=np.zeros(2)

######
# Here I start with a nested loop in order to explore all the points of the pes
for i in range(int(pointsy)):
    for j in range(int(pointsx)):
        if j>=1 and i>=1 and j<=(int(pointsx)-2) and i<=(int(pointsy)-2):
# This "if" statement lets us consider only the points that are NOT on the boundary of the pes,
# in fact, in those points the derivatives are not defined
            
            # Calculation of the Hessian
            H_xx=(energy[i,j+1]-2*energy[i,j]+energy[i,j-1])/(step_x**2)
            H_yy=(energy[i+1,j]-2*energy[i,j]+energy[i-1,j])/(step_y**2)
            H_xy=(energy[i+1,j+1]-energy[i+1,j-1]-energy[i-1,j+1]+energy[i-1,j-1])/(4*step_x*step_y)
            H_yx=H_xy
           
            Hessian=np.zeros((2,2))
            Hessian[0,0]=H_xx
            Hessian[1,1]=H_yy
            Hessian[1,0]=H_xy
            Hessian[0,1]=H_yx
            
            # Calculation of the matrices that contain the second derivatives of
            # of the dipole vector
            dipx_xx=(dipx[i,j+1]-2*dipx[i,j]+dipx[i,j-1])/(step_x**2)
            dipx_yy=(dipx[i+1,j]-2*dipx[i,j]+dipx[i-1,j])/(step_y**2)
            dipx_xy=(dipx[i+1,j+1]-dipx[i+1,j-1]-dipx[i-1,j+1]+dipx[i-1,j-1])/(4*step_x*step_y)
            dipx_yx=dipx_xy

            dipy_xx=(dipy[i,j+1]-2*dipy[i,j]+dipy[i,j-1])/(step_x**2)
            dipy_yy=(dipy[i+1,j]-2*dipy[i,j]+dipy[i-1,j])/(step_y**2)
            dipy_xy=(dipy[i+1,j+1]-dipy[i+1,j-1]-dipy[i-1,j+1]+dipy[i-1,j-1])/(4*step_x*step_y)
            dipy_yx=dipy_xy
            
            dipz_xx=(dipz[i,j+1]-2*dipz[i,j]+dipz[i,j-1])/(step_x**2)
            dipz_yy=(dipz[i+1,j]-2*dipz[i,j]+dipz[i-1,j])/(step_y**2)
            dipz_xy=(dipz[i+1,j+1]-dipz[i+1,j-1]-dipz[i-1,j+1]+dipz[i-1,j-1])/(4*step_x*step_y)
            dipz_yx=dipz_xy

            
            D1=np.zeros((2,2))
            D1[0,0]=dipx_xx
            D1[1,1]=dipx_yy
            D1[1,0]=dipx_xy
            D1[0,1]=dipx_yx
            
            D2=np.zeros((2,2))
            D2[0,0]=dipy_xx
            D2[1,1]=dipy_yy
            D2[1,0]=dipy_xy
            D2[0,1]=dipy_yx
            
            D3=np.zeros((2,2))
            D3[0,0]=dipz_xx
            D3[1,1]=dipz_yy
            D3[1,0]=dipz_xy
            D3[0,1]=D3[1,0]
           
            # Calculation of the gradient of the pes
            g_pes_x=(energy[i,j+1]-energy[i,j-1])/(2*step_x)
            g_pes_y=(energy[i+1,j]-energy[i-1,j])/(2*step_y)
            
            g_pes=np.zeros(2)
            g_pes[0]=g_pes_x
            g_pes[1]=g_pes_y
            
            
            
            # Calculation of the first derivatives of the dipole
            dipx_x=(dipx[i,j+1]-dipx[i,j-1])/(2*step_x)
            dipx_y=(dipx[i+1,j]-dipx[i-1,j])/(2*step_y)
            
            dipy_x=(dipy[i,j+1]-dipy[i,j-1])/(2*step_x)
            dipy_y=(dipy[i+1,j]-dipy[i-1,j])/(2*step_y)
            
            dipz_x=(dipz[i,j+1]-dipz[i,j-1])/(2*step_x)
            dipz_y=(dipz[i+1,j]-dipz[i-1,j])/(2*step_y)
            
            # I put the derivatives of the dipole in an array of the appropriate
            # shape
            grad_T_dipole=np.zeros((2,3))
            grad_T_dipole[0,0]=dipx_x
            grad_T_dipole[0,1]=dipy_x
            grad_T_dipole[0,2]=dipz_x
            grad_T_dipole[1,0]=dipx_y
            grad_T_dipole[1,1]=dipy_y
            grad_T_dipole[1,2]=dipz_y

            # Now I should calculate the P matrix (eq. 8 ref. work), we found out
            # that is numerically more stable to compute it using the SciPy library
            #  than computing it explicitly according to eq. 8
            P_2=sp.linalg.pinv(grad_T_dipole)
            print('Moore-Penrose')
            print(P_2)
            # With this "if" statement I check two things: 1) The existance of 
            # the P matrix and 2) that the gradient of the unpertubed pes
            # is non-zero
            test_det=np.dot(grad_T_dipole,P_2)
            print(np.linalg.det(test_det))
            
            if np.abs(np.linalg.det(test_det))>1e-5 and np.linalg.norm(g_pes)>1e-5:
                # Calculation of e*E=P*g and storage in "fields" array
                Pg=np.dot(P_2,g_pes)
                fields=fields=np.vstack((fields,Pg))
                
                # Calculation and storage of the perturbed hessian (eq 12 ref. work)
                perturbed_hessian=Hessian-D1*Pg[0]-D2*Pg[1]-D3*Pg[2]
                H_e=np.vstack((H_e,perturbed_hessian[None]))
                
                # Calculation and storage of s and of the norm of the "s" vector (eq 22 ref. work)
                norm_s=np.append(norm_s,np.linalg.norm(np.dot(perturbed_hessian,g_pes)))
                tmp_s=np.dot(perturbed_hessian,g_pes)
                s=np.vstack((s,tmp_s))
                
                # storage of the indexes of the foor loop
                # needed to retrieve the geometry of the OBBP
                tmp=np.array([i,j])
                
                indexes=np.vstack((indexes,tmp))
                
                # Calculation of the perturbed gradient for each point of the grid
                tmp_perturbed_gradient=g_pes-np.dot(grad_T_dipole,Pg)
                perturbed_gradient=np.vstack((perturbed_gradient,tmp_perturbed_gradient))
                
                # Storage of the gradient of the original pes
                gradient_pes=np.vstack((gradient_pes,g_pes))

# Cleaning of the storage arrays 
gradient_pes = np.delete(gradient_pes, (0), axis=0)
perturbed_gradient = np.delete(perturbed_gradient, (0), axis=0)
fields = np.delete(fields, (0), axis=0)   
s = np.delete(s, (0), axis=0)    
H_e = np.delete(H_e, (0), axis=0)            
indexes = np.delete(indexes, (0), axis=0)                 
indexes=indexes.astype(int)                    

# A couple of words about the storage of H_e. I choose to store each H_e matrix
# in a 3D array in which every "layer" contains the matrix associated to a point 
# of the grid (if it passed the tests regarding the existance of P and the value
# of the gradient). Remember the indexing for 3D arrays: it is not, naively,
# i,j,k= rows, columns,layers but i,j,k=layers, rows, columns



print('---- RESULTS ----')
print('')
# For each point of the pes I check if the unperturbed gradient is an eigenvector
# with eigenvalue zero of the perturbed hessian. To do so I normalize the gradient 
# and multiply the pertubed hessian and the normalized gradient. If it is an eigenvector
# with eigenvalue zero the norm of the vector resulting from the product should be zero.
# I check if this condition is true for every point of the grid. 
# Numerically only the points of the grid that give a norm of the vector "H_e*gradient_normalized"
# less than a treshold will be further examined. The value of the aforementioned norm will be 
# printed as "EIGENVECTOR TEST" only if other criteria will be met.
for i in range(len(norm_s)):
               normalized_gradient=(gradient_pes[i,:])/np.linalg.norm(gradient_pes[i,:])
               a=np.dot(H_e[i,:,:],normalized_gradient)
               # print(np.linalg.norm(a))
               if np.linalg.norm(a)<2e-2:
# Now that the gradient is an eigenvector with eigenvalue zero of H_e, I need to check
# the geometry associated to the point. To do so it is necessary to do some index gymnastics.
# I stored in "indexes" the indexes of the points of the grid that passed the previous test and
# now I check to which geometry they correspond looking at the energy.
                   tmp_index=(indexes[i])
                   for j in range(len(energy_column)):
                       if energy[tmp_index[0],tmp_index[1]]==energy_column[j]:
                           print('geometry_index=',j)
                           print('P-C=',data[j,0])
                           print('C-O=',data[j,1])
# I found the point that satisfies the requirements for being an OBBP, but I still have to check
# if it belongs to the minimum energy path, that is I need to check if the energy of the 
# point is lower than the one of the transition state. If this is true I print the aforementioned 
# "EIGENVECTOR TEST", the value of the norm of the "s" vector, the perturbed gradient,
# the corresponding field and, as a final test, the determinant of the H_e matrix.
                           print('IS THE ENERGY OF THIS POINT LOWER THAN THE ONE OF THE TS?')
                           if float(energy_column[j])<energy_ts:
                                print('YES')
                                print('ENERGY TS - ENERGY THIS POINT? (kcal/mol)')
                                print((energy_ts-energy_column[j])*627.503)
                                print('ENERGY OP1 - ENERGY THIS POINT? (kcal/mol)')
                                print((energy_op1-energy_column[j])*627.503)  
                                print('EIGENVECTOR TEST (norm of H_e*g_normalized)')
                                print(np.linalg.norm(a))
                                b=np.dot(np.transpose(gradient_pes[i,:]),H_e[i,:,:])
                                print('g_transpose*H_e*g')
                                print(np.dot(b,gradient_pes[i,:]))
                                print('VALUE NORM s=H_e*g')
                                print(norm_s[i])
                                print('IS DET(H_e)=0 IN THIS POINT?')
                                print(np.linalg.det(H_e[i,:,:]))
                               
                                print('CORRESPONDING PERTURBED GRADIENT')
                                print(perturbed_gradient[i,:])
                                print('CORRESPONDING UNPERTURBED GRADIENT')
                                print(gradient_pes[i,:])
                                print('CORRESPONDING FIELD (a.u.)')
                                print(fields[i,:])
                                print('NORM OF THE FIELD (a.u.)')
                                norm_field=np.linalg.norm(fields[i,:])
                                print(norm_field)
                                print('NORMALIZED FIELD (a.u.)')
                                print(fields[i,:]/(norm_field))
                                print('CORRESPONDING FIELD (V/m)')
                                print(fields[i,:]* 5.14220826*10**11)
                                print(tmp_index)
                           else:
                                print('NO')
                                print('energy this point =',energy_column[j])
                                print('energy TS =',energy_ts)
                           
                           print('#######################')
                           print('#######################')
                           print('')
                           print('')




