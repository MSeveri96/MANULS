#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 08:48:58 2022

@author: Marco Severi
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   :  14}

matplotlib.rc('font', **font)

######################################################################
##################### USER DEFINED SECTION  ##########################
######################################################################

# This is the only section that needs to be modified, the rest is up to the code

# Loading of the files containing the data, please see the documentation

data=np.loadtxt('/home/marco/huisgen/tzvp/from_reacts/scan_tzvp_from_reacts_grid.txt')
polar=np.loadtxt('/home/marco/huisgen/tzvp/from_reacts/polarizability_scan_tzvp_from_reacts.txt')


# "pointsx" and "pointsy" = number of points scanned for each variable
pointsx=int(25)
pointsy=int(25)


# "step_x" and "step_y" = step size used in the scan 
step_x=0.025
step_y=0.025

# number of directions scanned in the search of the optimal field
# usually the default gives a satisfactory time/precision ratio
n_directions=100000

#geometry of the reactants
geometry_reactant_x=3.5
geometry_reactant_y=3.5


    

######################################################################
#####################      WARNING!!!     ##########################
######################################################################
# Given the construction of the code it is very important to put the data in a 
# precise order. Example: Suppose that you want to find the optimal field that
# breaks two bonds. The scan must be performed starting at short bond lengths 
# and finishing at large bond lengths. If the data are not supplied in this
# format the location of the optimal BBP and the derivatives will be 
# incorrect

######################################################################
##################### END OF USER DEFINED SECTION  ###################
######################################################################




########################################################################
########################################################################
########################################################################


print('THE CODE IS WORKING...')

f_tol=5e-5
m_tol=0.999
g_pert_tol=5e-5
energy_tol=-0.0048

indexes=np.zeros(2)
gradient_pes=np.zeros(2)
vector_norm_gradient_extr=np.array([])
H=np.zeros((1,2,2))

optimal_f=np.array([])
perturbed_gradient_optimal_fields=np.zeros(2)
original_gradient_optimal_field=np.zeros(2)
optimal_fields=np.zeros(3)
norm_optimal_fields=np.array([])



                    
            



# Starting from the data written in columns I rearrange them in matrices, in 
# order to implement easily the derivatives. Also in this case the ordering
# "Fortran style" is very important. At the end the columns reprensent
#  the x-axis, while the rows represent the y-axis


energy=np.reshape(data[:,2],(pointsy,pointsx),order='F')
dipx=np.reshape(data[:,3],(pointsy,pointsx),order='F')
dipy=np.reshape(data[:,4],(pointsy,pointsx),order='F')
dipz=np.reshape(data[:,5],(pointsy,pointsx),order='F')

energy_column=data[:,2]


# Reading and rearranging the polarizability


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


# Here the code finds the point on the supplied grid that is nearest to the 
#geometry of the reactants. 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


if geometry_reactant_x!=-1000:
    geom_react_on_grid_x=find_nearest(data[:,0], geometry_reactant_x)
    geom_react_on_grid_y=find_nearest(data[:,1], geometry_reactant_y)
    
    for row in range((pointsy)):
        for column in range((pointsx)):
            for idx_data in range(len(data[:,0])):
                if energy[row,column]==data[idx_data,2] and data[idx_data,0]==geom_react_on_grid_x and data[idx_data,1]==geom_react_on_grid_y:
                    row_react=int(row)
                    column_react=int(column)

file_optimal= open("external_electric_fields.txt","w+")
file_bbps= open("coordinates_bbps.txt","w+")


# The code finds the optimal BBP. First it computtes the gradient and the 
# hessian for each point of the grid, then it check the gradient extremal 
# condition. Between the points in which the condition is satisfied it
# pick as optimal BBP the point in which the gradien is maximum.

for i in range(int(pointsy)):
        for j in range(int(pointsx)):
            if j>=1 and i>=1 and j<=(int(pointsx)-2) and i<=(int(pointsy)-2):
                g_pes_x=(energy[i,j+1]-energy[i,j-1])/(2*step_x)
                g_pes_y=(energy[i+1,j]-energy[i-1,j])/(2*step_y)
                
                g_pes=np.zeros(2)
                g_pes[0]=g_pes_x
                g_pes[1]=g_pes_y
                tmp=np.array([i,j])
                    
                indexes=np.vstack((indexes,tmp))
                gradient_pes=np.vstack((gradient_pes,g_pes))
                
                H_xx=(energy[i,j+1]-2*energy[i,j]+energy[i,j-1])/(step_x**2)
                H_yy=(energy[i+1,j]-2*energy[i,j]+energy[i-1,j])/(step_y**2)
                H_xy=(energy[i+1,j+1]-energy[i+1,j-1]-energy[i-1,j+1]+energy[i-1,j-1])/(4*step_x*step_y)
                H_yx=H_xy
           
                Hessian=np.zeros((2,2))
                Hessian[0,0]=H_xx
                Hessian[1,1]=H_yy
                Hessian[1,0]=H_xy
                Hessian[0,1]=H_yx
                H=np.vstack((H,Hessian[None]))
                




indexes = np.delete(indexes, (0), axis=0)
gradient_pes = np.delete(gradient_pes, (0), axis=0)
H = np.delete(H, (0), axis=0)



for idx_ge in range(len(indexes)):
    if gradient_pes[idx_ge,0]>0 or gradient_pes[idx_ge,1]>0:
        norm_sq_g_pes=np.matmul(gradient_pes[idx_ge,:],gradient_pes[idx_ge,:])
        I=np.eye(2)
        tmp_proj=np.outer(gradient_pes[idx_ge,:],gradient_pes[idx_ge,:])
        proj=I-(tmp_proj/norm_sq_g_pes)
        tmp_grad_ext=np.matmul(proj,H[idx_ge,:,:])
        gradient_extremal=np.matmul(tmp_grad_ext,gradient_pes[idx_ge,:])
        vector_norm_gradient_extr=np.append(vector_norm_gradient_extr,np.linalg.norm(gradient_extremal))


gradient_maxima = np.argsort(vector_norm_gradient_extr)[::1][:5]
array_for_max_gradient=np.array([])


for l in range(len(gradient_maxima)):
    norm_gradient=np.linalg.norm(gradient_pes[gradient_maxima[l],:])
    array_for_max_gradient=np.append(array_for_max_gradient,norm_gradient)
    obbp_idx=indexes[gradient_maxima[l]]
    obbp_x=int(obbp_idx[1])
    obbp_y=int(obbp_idx[0])
    for c in range(len(energy_column)):
        if energy[obbp_y,obbp_x]==energy_column[c]:
            print('Geometry of the BBP',file=file_bbps)
            print('variable 1 =',data[c,0],file=file_bbps)
            print('variable 2 =',data[c,1],file=file_bbps)
            print('index BBP =',indexes[gradient_maxima[l]],file=file_bbps)
            print('norm gradient in this point=',norm_gradient,file=file_bbps)
            print('gradient in this point',gradient_pes[gradient_maxima[l]],file=file_bbps)
            print('gradient extremal condition',vector_norm_gradient_extr[gradient_maxima[l]],file=file_bbps)
            print('#########',file=file_bbps)
            print('',file=file_bbps)
            print('',file=file_bbps)




obbp_idx=indexes[gradient_maxima[np.argmax(array_for_max_gradient)]]
obbp_x=int(obbp_idx[1])
obbp_y=int(obbp_idx[0])



for c in range(len(energy_column)):
    if energy[obbp_y,obbp_x]==energy_column[c]:
        print('Geometry of the optimal BBP',file=file_optimal)
        print('variable 1 =',data[c,0],file=file_optimal)
        print('variable 2 =',data[c,1],file=file_optimal)
        # print('Geometry of the optimal BBP')
        # print('variable 1 =',data[c,0])
        # print('variable 2 =',data[c,1])


print('gradient extremal condition',vector_norm_gradient_extr[gradient_maxima[np.argmax(array_for_max_gradient)]],file=file_optimal)
print('gradient extremal condition',vector_norm_gradient_extr[gradient_maxima[np.argmax(array_for_max_gradient)]])
print('norm original gradient',np.amax(array_for_max_gradient),file=file_optimal)
print('--------------------------------------------------------',file=file_optimal)
print('--------------------------------------------------------',file=file_optimal)


# Now the optimal BBP is defined. The code computes the ingredients for
# finding the optimal field. It starts computing the hessian and
# the first and second derivatives of the dipole moment components

dipx_x=(dipx[obbp_y,obbp_x+1]-dipx[obbp_y,obbp_x-1])/(2*step_x)
dipx_y=(dipx[obbp_y+1,obbp_x]-dipx[obbp_y-1,obbp_x])/(2*step_y)
        
dipy_x=(dipy[obbp_y,obbp_x+1]-dipy[obbp_y,obbp_x-1])/(2*step_x)
dipy_y=(dipy[obbp_y+1,obbp_x]-dipy[obbp_y-1,obbp_x])/(2*step_y)
        
dipz_x=(dipz[obbp_y,obbp_x+1]-dipz[obbp_y,obbp_x-1])/(2*step_x)
dipz_y=(dipz[obbp_y+1,obbp_x]-dipz[obbp_y-1,obbp_x])/(2*step_y)
        
grad_T_dipole=np.zeros((2,3))
grad_T_dipole[0,0]=dipx_x
grad_T_dipole[0,1]=dipy_x
grad_T_dipole[0,2]=dipz_x
grad_T_dipole[1,0]=dipx_y
grad_T_dipole[1,1]=dipy_y
grad_T_dipole[1,2]=dipz_y



g_pes_x_obbp=(energy[obbp_y,obbp_x+1]-energy[obbp_y,obbp_x-1])/(2*step_x)
g_pes_y_obbp=(energy[obbp_y+1,obbp_x]-energy[obbp_y-1,obbp_x])/(2*step_y)
            
g_pes_obbp=np.zeros(2)
g_pes_obbp[0]=g_pes_x_obbp
g_pes_obbp[1]=g_pes_y_obbp



H_obbp_xx=(energy[obbp_y,obbp_x+1]-2*energy[obbp_y,obbp_x]+energy[obbp_y,obbp_x-1])/(step_x**2)
H_obbp_yy=(energy[obbp_y+1,obbp_x]-2*energy[obbp_y,obbp_x]+energy[obbp_y-1,obbp_x])/(step_y**2)
H_obbp_xy=(energy[obbp_y+1,obbp_x+1]-energy[obbp_y+1,obbp_x-1]-energy[obbp_y-1,obbp_x+1]+energy[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
H_obbp_yx=H_obbp_xy

Hessian_obbp=np.zeros((2,2))
Hessian_obbp[0,0]=H_obbp_xx
Hessian_obbp[1,1]=H_obbp_yy
Hessian_obbp[1,0]=H_obbp_xy
Hessian_obbp[0,1]=H_obbp_yx

dipx_xx=(dipx[obbp_y,obbp_x+1]-2*dipx[obbp_y,obbp_x]+dipx[obbp_y,obbp_x-1])/(step_x**2)
dipx_yy=(dipx[obbp_y+1,obbp_x]-2*dipx[obbp_y,obbp_x]+dipx[obbp_y-1,obbp_x])/(step_y**2)
dipx_xy=(dipx[obbp_y+1,obbp_x+1]-dipx[obbp_y+1,obbp_x-1]-dipx[obbp_y-1,obbp_x+1]+dipx[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
dipx_yx=dipx_xy

dipy_xx=(dipy[obbp_y,obbp_x+1]-2*dipy[obbp_y,obbp_x]+dipy[obbp_y,obbp_x-1])/(step_x**2)
dipy_yy=(dipy[obbp_y+1,obbp_x]-2*dipy[obbp_y,obbp_x]+dipy[obbp_y-1,obbp_x])/(step_y**2)
dipy_xy=(dipy[obbp_y+1,obbp_x+1]-dipy[obbp_y+1,obbp_x-1]-dipy[obbp_y-1,obbp_x+1]+dipy[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
dipy_yx=dipy_xy

dipz_xx=(dipz[obbp_y,obbp_x+1]-2*dipz[obbp_y,obbp_x]+dipz[obbp_y,obbp_x-1])/(step_x**2)
dipz_yy=(dipz[obbp_y+1,obbp_x]-2*dipz[obbp_y,obbp_x]+dipz[obbp_y-1,obbp_x])/(step_y**2)
dipz_xy=(dipz[obbp_y+1,obbp_x+1]-dipz[obbp_y+1,obbp_x-1]-dipz[obbp_y-1,obbp_x+1]+dipz[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
dipz_yx=dipz_xy


# The code generates "n_directions" unit vectors using the Fibonacci
# sphere algorithm


def fibonacci_sphere(samples):

    points = []
    phi = math.pi * (3. - math.sqrt(5.)) 

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  
        radius = math.sqrt(1 - y * y)  

        theta = phi * i  

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points

points=fibonacci_sphere(n_directions)

points=np.array(points)


# For each generated field direction the code computes the 
# remaining "ingredients" for the analysis.

for g in range(len(points)):
    en=np.array([points[g,0],points[g,1],points[g,2]])
    
    
    

    fe0=np.zeros((pointsy,pointsx))
    fe1=np.zeros((pointsy,pointsx))
    fe2=np.zeros((pointsy,pointsx))
    for h in range(pointsy):
        for f in range(pointsx):
            fe0[h,f]=en[0]*a00[h,f]+en[1]*a01[h,f]+en[2]*a02[h,f]
            fe1[h,f]=en[0]*a10[h,f]+en[1]*a11[h,f]+en[2]*a12[h,f]
            fe2[h,f]=en[0]*a20[h,f]+en[1]*a21[h,f]+en[2]*a22[h,f]
            
    fe0_x=(fe0[obbp_y,obbp_x+1]-fe0[obbp_y,obbp_x-1])/(2*step_x)
    fe0_y=(fe0[obbp_y+1,obbp_x]-fe0[obbp_y-1,obbp_x])/(2*step_y)
            
    fe1_x=(fe1[obbp_y,obbp_x+1]-fe1[obbp_y,obbp_x-1])/(2*step_x)
    fe1_y=(fe1[obbp_y+1,obbp_x]-fe1[obbp_y-1,obbp_x])/(2*step_y)
            
    fe2_x=(fe2[obbp_y,obbp_x+1]-fe2[obbp_y,obbp_x-1])/(2*step_x)
    fe2_y=(fe2[obbp_y+1,obbp_x]-fe2[obbp_y-1,obbp_x])/(2*step_y)
            
    
    grad_T_dipole=np.zeros((2,3))
    grad_T_dipole[0,0]=dipx_x
    grad_T_dipole[0,1]=dipy_x
    grad_T_dipole[0,2]=dipz_x
    grad_T_dipole[1,0]=dipx_y
    grad_T_dipole[1,1]=dipy_y
    grad_T_dipole[1,2]=dipz_y

    grad_T_fe=np.zeros((2,3))
    grad_T_fe[0,0]=fe0_x
    grad_T_fe[0,1]=fe1_x
    grad_T_fe[0,2]=fe2_x
    grad_T_fe[1,0]=fe0_y
    grad_T_fe[1,1]=fe1_y
    grad_T_fe[1,2]=fe2_y
    
    fe0_xx=(fe0[obbp_y,obbp_x+1]-2*fe0[obbp_y,obbp_x]+fe0[obbp_y,obbp_x-1])/(step_x**2)
    fe0_yy=(fe0[obbp_y+1,obbp_x]-2*fe0[obbp_y,obbp_x]+fe0[obbp_y-1,obbp_x])/(step_y**2)
    fe0_xy=(fe0[obbp_y+1,obbp_x+1]-fe0[obbp_y+1,obbp_x-1]-fe0[obbp_y-1,obbp_x+1]+fe0[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
    fe0_yx=fe0_xy
    
    second_der_f0=np.zeros((2,2))
    second_der_f0[0,0]=fe0_xx
    second_der_f0[1,1]=fe0_yy
    second_der_f0[1,0]=fe0_xy
    second_der_f0[0,1]=fe0_yx
    
    fe1_xx=(fe1[obbp_y,obbp_x+1]-2*fe1[obbp_y,obbp_x]+fe1[obbp_y,obbp_x-1])/(step_x**2)
    fe1_yy=(fe1[obbp_y+1,obbp_x]-2*fe1[obbp_y,obbp_x]+fe1[obbp_y-1,obbp_x])/(step_y**2)
    fe1_xy=(fe1[obbp_y+1,obbp_x+1]-fe1[obbp_y+1,obbp_x-1]-fe1[obbp_y-1,obbp_x+1]+fe1[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
    fe1_yx=fe1_xy
    
    second_der_f1=np.zeros((2,2))
    second_der_f1[0,0]=fe1_xx
    second_der_f1[1,1]=fe1_yy
    second_der_f1[1,0]=fe1_xy
    second_der_f1[0,1]=fe1_yx
    
    fe2_xx=(fe2[obbp_y,obbp_x+1]-2*fe2[obbp_y,obbp_x]+fe2[obbp_y,obbp_x-1])/(step_x**2)
    fe2_yy=(fe2[obbp_y+1,obbp_x]-2*fe2[obbp_y,obbp_x]+fe2[obbp_y-1,obbp_x])/(step_y**2)
    fe2_xy=(fe2[obbp_y+1,obbp_x+1]-fe2[obbp_y+1,obbp_x-1]-fe2[obbp_y-1,obbp_x+1]+fe2[obbp_y-1,obbp_x-1])/(4*step_x*step_y)
    fe2_yx=fe2_xy
    
    second_der_f2=np.zeros((2,2))
    second_der_f2[0,0]=fe2_xx
    second_der_f2[1,1]=fe2_yy
    second_der_f2[1,0]=fe2_xy
    second_der_f2[0,1]=fe2_yx
    
    h_en=np.dot(grad_T_dipole,en)
    
    p_en=np.dot(grad_T_fe,en)
        
# The code computes the amplitude of the electric field in the given
# direction
    coeff=np.zeros(5)
    coeff[4]=(1/4)*(np.dot(p_en,p_en))
    coeff[3]=(np.dot(h_en,p_en))
    coeff[2]=(np.dot(h_en,h_en))
    coeff[1]=0.0000000
    coeff[0]=-(np.dot(g_pes_obbp,g_pes_obbp))
    
    
    amplitudes=np.polynomial.polynomial.polyroots(coeff)
    amplitudes=np.array(amplitudes)
    
    real_valued = amplitudes.real[abs(amplitudes.imag)<1e-8]
    

    if not  real_valued.size == 0:
        idx_min_sol=np.where(real_valued > 0, real_valued, np.inf).argmin()
        
        E=real_valued[idx_min_sol] 
        
    else:
        E=100
    
    e=E*en

    Mx=np.zeros((2,2))
    Mx[0,0]=dipx_xx+(1/2)*fe0_xx*(E**2)
    Mx[1,1]=dipx_yy+(1/2)*fe0_yy*(E**2)
    Mx[1,0]=dipx_xy+(1/2)*fe0_xy*(E**2)
    Mx[0,1]=dipx_yx+(1/2)*fe0_yx*(E**2)
    
    
    My=np.zeros((2,2))
    My[0,0]=dipy_xx+(1/2)*fe1_xx*(E**2)
    My[1,1]=dipy_yy+(1/2)*fe1_yy*(E**2)
    My[1,0]=dipy_xy+(1/2)*fe1_xy*(E**2)
    My[0,1]=dipy_yx+(1/2)*fe1_yx*(E**2)
    
    Mz=np.zeros((2,2))
    Mz[0,0]=dipz_xx+(1/2)*fe2_xx*(E**2)
    Mz[1,1]=dipz_yy+(1/2)*fe2_yy*(E**2)
    Mz[1,0]=dipz_xy+(1/2)*fe2_xy*(E**2)
    Mz[0,1]=dipz_yx+(1/2)*fe2_yx*(E**2)
    
    matrix_of_gradients=grad_T_dipole+(1/2)*E*grad_T_fe
    
    D=np.matmul(grad_T_dipole.transpose(),grad_T_dipole)
    normalized_gradient_obbp=g_pes_obbp/np.linalg.norm(g_pes_obbp)
    m=np.matmul(matrix_of_gradients.transpose(),normalized_gradient_obbp)
    D_bar=D-np.outer(m,m)
    
    norm_sq_g_pes_obbp=np.matmul(g_pes_obbp.transpose(),g_pes_obbp)
    I=np.eye(2)
    tmp_proj_obbp=np.outer(g_pes_obbp,g_pes_obbp)
    proj_obbp=I-(tmp_proj_obbp/norm_sq_g_pes_obbp)
    
    tmp_v_vec=np.matmul(proj_obbp,Hessian_obbp)
    v_vec=np.matmul(tmp_v_vec,g_pes_obbp)
    
    
    tmp_test=np.matmul(g_pes_obbp.transpose(),Hessian_obbp)
    v_italic=np.matmul(tmp_test,g_pes_obbp)
    
    
    
    tx=np.dot(Mx,g_pes_obbp)
    ty=np.dot(My,g_pes_obbp)
    tz=np.dot(Mz,g_pes_obbp)
    
    wx=np.dot(g_pes_obbp,tx)
    wy=np.dot(g_pes_obbp,ty)
    wz=np.dot(g_pes_obbp,tz)
    w_vec=np.array([wx,wy,wz])
    
    
    
    
    Z2=np.zeros((3,3))
    
    Z2[0,0]=np.matmul(tx,tx)
    Z2[0,1]=np.matmul(tx,ty)
    Z2[0,2]=np.matmul(tx,tz)
    Z2[1,0]=np.matmul(ty,tx)
    Z2[1,1]=np.matmul(ty,ty)
    Z2[1,2]=np.matmul(ty,tz)
    Z2[2,0]=np.matmul(tz,tx)
    Z2[2,1]=np.matmul(tz,ty)
    Z2[2,2]=np.matmul(tz,tz)
    
# All the ingredients are ready. The code computes the f function

    tmp_eq1=np.dot(e,D_bar)
    tmp_eq_3_1=np.matmul(e,Z2)
    tmp_eq_3_2=np.matmul(v_vec,v_vec)+v_italic**2
    
    b=np.array([0,v_italic,tmp_eq_3_2])
    
    tmp_B=np.row_stack((tmp_eq1,w_vec))
    B=np.row_stack((tmp_B,tmp_eq_3_1))
    
    tmp_f_1=np.matmul(e.transpose(),B)
    tmp_f_2=np.matmul(b.transpose(),B)
    eq2=0
    
    f=(1/2)*np.matmul(tmp_f_1,e)-np.matmul(tmp_f_2,e)+(1/2)*np.matmul(b.transpose(),b)
    
    energy_perturbed=np.zeros((pointsy,pointsx))

#  The code computes the energy perturbed by the fied and 
#  performs the test (see documentation for further details)

    for idx_p_y in range(int(pointsy)):
        for idx_p_x in range(int(pointsx)):
            energy_perturbed[idx_p_y,idx_p_x]=energy[idx_p_y,idx_p_x]-(dipx[idx_p_y,idx_p_x]*e[0]+dipy[idx_p_y,idx_p_x]*e[1]+dipz[idx_p_y,idx_p_x]*e[2])-(1/2)*e[0]*(e[0]*a00[idx_p_y,idx_p_x]+e[1]*a01[idx_p_y,idx_p_x]+e[2]*a02[idx_p_y,idx_p_x])-(1/2)*e[1]*(e[0]*a10[idx_p_y,idx_p_x]+e[1]*a11[idx_p_y,idx_p_x]+e[2]*a12[idx_p_y,idx_p_x])-(1/2)*e[2]*(e[0]*a20[idx_p_y,idx_p_x]+e[1]*a21[idx_p_y,idx_p_x]+e[2]*a22[idx_p_y,idx_p_x])
            
    
    if geometry_reactant_x!=-1000:
        
        if np.abs(f)<f_tol and np.abs(np.dot(m,e)/np.linalg.norm(g_pes_obbp))>m_tol and np.linalg.norm(g_pes_obbp-np.dot(matrix_of_gradients,e))<g_pert_tol and energy_perturbed[row_react,column_react]-energy_perturbed[obbp_y,obbp_x]>energy_tol:
            # print(f)
            # print(np.abs(np.dot(m,e)/np.linalg.norm(g_pes_obbp)))
            # print(np.linalg.norm(g_pes_obbp-np.dot(matrix_of_gradients,e)))
            # print(np.linalg.norm(g_pes_obbp))
            # print('field',e)
            # print(np.linalg.norm(e))
            # print(e/np.linalg.norm(e))
            # print('pert function', energy_perturbed[obbp_y,obbp_x]-energy[obbp_y,obbp_x])
            
            # print('----')   
            
            optimal_f=np.append(optimal_f,f)
            perturbed_gradient_optimal_fields=np.vstack((perturbed_gradient_optimal_fields,g_pes_obbp-np.dot(matrix_of_gradients,e)))
            original_gradient_optimal_field=np.vstack((original_gradient_optimal_field,g_pes_obbp))
            optimal_fields=np.vstack((optimal_fields,e))
            norm_optimal_fields=np.append(norm_optimal_fields,np.linalg.norm(e))
    else:
        if np.abs(f)<f_tol and np.abs(np.dot(m,e)/np.linalg.norm(g_pes_obbp))>m_tol and np.linalg.norm(g_pes_obbp-np.dot(matrix_of_gradients,e))<g_pert_tol:
            # print(f)
            # print(np.abs(np.dot(m,e)/np.linalg.norm(g_pes_obbp)))
            # print(np.linalg.norm(g_pes_obbp-np.dot(matrix_of_gradients,e)))
            # print(np.linalg.norm(g_pes_obbp))
            # print('field',e)
            # print(np.linalg.norm(e))
            # print(e/np.linalg.norm(e))
            # print('----')   
            
            optimal_f=np.append(optimal_f,f)
            perturbed_gradient_optimal_fields=np.vstack((perturbed_gradient_optimal_fields,g_pes_obbp-np.dot(matrix_of_gradients,e)))
            original_gradient_optimal_field=np.vstack((original_gradient_optimal_field,g_pes_obbp))
            optimal_fields=np.vstack((optimal_fields,e))
            norm_optimal_fields=np.append(norm_optimal_fields,np.linalg.norm(e))

perturbed_gradient_optimal_fields = np.delete(perturbed_gradient_optimal_fields, (0), axis=0)
original_gradient_optimal_field = np.delete(original_gradient_optimal_field, (0), axis=0)
optimal_fields = np.delete(optimal_fields, (0), axis=0)



# The code did the things above for all the "n_directions" and stored all
# the fields that met the requirements. Now it picks the field with the 
# lowest amplitude as optimal.


idx_optimal_field=np.argsort(norm_optimal_fields,kind='stable')

for z in range(len(idx_optimal_field)):
    if z==0:
        print('Results:')
        print('')
        print('')
        print('Optimal electric field (a.u.)',file=file_optimal)
        print(optimal_fields[idx_optimal_field[z],:],file=file_optimal)
        print('Amplitude of the optimal electric field (a.u.)',file=file_optimal)
        
        print(norm_optimal_fields[idx_optimal_field[z]],file=file_optimal)
        print('Direction of the optimal electric field',file=file_optimal)
        print(optimal_fields[idx_optimal_field[z],:]/norm_optimal_fields[idx_optimal_field[z]],file=file_optimal)
        
        print('Optimal electric field (V/m)',file=file_optimal)
        print(optimal_fields[idx_optimal_field[z],:]* 5.14220826*10**11,file=file_optimal)
        
        print('Amplitude of the optimal electric field (V/m)',file=file_optimal)
        print('{:.4e}'.format(norm_optimal_fields[idx_optimal_field[z]]* 5.14220826*10**11),file=file_optimal)
        
        print('Original Gradient',file=file_optimal)
        print(original_gradient_optimal_field[idx_optimal_field[z],:],file=file_optimal)
        print('Perturbed Gradient',file=file_optimal)
        print(perturbed_gradient_optimal_fields[idx_optimal_field[z],:],file=file_optimal)
        print('f Function',file=file_optimal)
        print(optimal_f[idx_optimal_field[z]],file=file_optimal)
        print('------------------------',file=file_optimal)
        print('',file=file_optimal)
        print('',file=file_optimal)
    if z>0:
        print('Sub-optimal electric field (a.u.)',file=file_optimal)
        print(optimal_fields[idx_optimal_field[z],:],file=file_optimal)
        print('Amplitude of the sub-optimal electric field (a.u.)',file=file_optimal)
        
        print(norm_optimal_fields[idx_optimal_field[z]],file=file_optimal)
        print('Direction of the sub-optimal electric field',file=file_optimal)
        print(optimal_fields[idx_optimal_field[z],:]/norm_optimal_fields[idx_optimal_field[z]],file=file_optimal)
        
        print('Sub-optimal electric field (V/m)',file=file_optimal)
        print(optimal_fields[idx_optimal_field[z],:]* 5.14220826*10**11,file=file_optimal)
        
        print('Amplitude of the sub-optimal electric field (V/m)',file=file_optimal)
        print('{:.4e}'.format(norm_optimal_fields[idx_optimal_field[z]]* 5.14220826*10**11),file=file_optimal)
        
        print('Original Gradient',file=file_optimal)
        print(original_gradient_optimal_field[idx_optimal_field[z],:],file=file_optimal)
        print('Perturbed Gradient',file=file_optimal)
        print(perturbed_gradient_optimal_fields[idx_optimal_field[z],:],file=file_optimal)
        print('f Function',file=file_optimal)
        print(optimal_f[idx_optimal_field[z]],file=file_optimal)
        
        print('',file=file_optimal)
        print('',file=file_optimal)





            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            






