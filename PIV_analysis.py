# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 11:09:42 2021

@author: robin
"""

import numpy as np
import os
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from scipy import signal
import math
import pandas as pd

rows=100
columns=100
zones=500

U=np.zeros((rows,columns,zones))
V=np.zeros((rows,columns,zones))
Y=np.zeros((rows,columns,zones))
Uflc=np.zeros((rows,columns,zones))
Uflcp=np.zeros((rows,columns,zones))
Uflcn=np.zeros((rows,columns,zones))
Vflcp=np.zeros((rows,columns,zones))
Vflcn=np.zeros((rows,columns,zones))
Vflc=np.zeros((rows,columns,zones))
Umean=np.zeros((rows,columns))
Vmean=np.zeros((rows,columns))
Ymean=np.zeros((rows,columns))
Urms=np.zeros((rows,columns))
Urmsp=np.zeros((rows,columns))
Vrmsp=np.zeros((rows,columns))
Vrms=np.zeros((rows,columns))
Urmsn=np.zeros((rows,columns))
Vrmsn=np.zeros((rows,columns))
loc=int(input('Enter the location at which you want to plot'))
file_path= "E:\OneDrive - Indian Institute of Science\Codes\PIV analysis"
data=np.loadtxt('X325_YUV_clean')
m=0

print('Writing velocity')
for k in range(0,zones):

    for i in range(0,columns):
        for j in range(0,rows):
            Y[i,j,k]=data[m,0]
            U[i,j,k]=data[m,1]
            V[i,j,k]=data[m,2]
            m=m+1

for i in range(0,columns):
    for j in range(0,rows):
        Ymean[i,j]=np.mean(Y[i,j,:])
        Umean[i,j]=np.mean(U[i,j,:])
        Vmean[i,j]=np.mean(V[i,j,:])

print('Calculating fluctuations')
for k in range(0,zones):

    for i in range(0,columns):
        for j in range(0,rows):
            Uflc[i,j,k]=U[i,j,k]-Umean[i,j]
            Vflc[i,j,k]=V[i,j,k]-Vmean[i,j]

for i in range(0,columns):
    for j in range(0,rows):
        Urms[i,j]=np.std(Uflc[i,j,:])
        Vrms[i,j]=np.std(Vflc[i,j,:])

print('Separating +ve and -ve')
for k in range(0,zones):

    for i in range(0,columns):
        for j in range(0,rows):
            if(Uflc[i,j,k]>0):
               Uflcp[i,j,k]=Uflc[i,j,k]
            elif(Uflc[i,j,k]<0):
               Uflcn[i,j,k]=Uflc[i,j,k]

for i in range(0,columns):
    for j in range(0,rows):
        tempp=Uflcp[i,j,:]
        temp1=tempp[np.where(tempp!=0)]
        Urmsp[i,j]=np.sqrt(np.mean(temp1**2))
        tempn=Uflcn[i,j,:]
        temp2=tempn[np.where(tempn!=0)]
        Urmsn[i,j]=np.sqrt(np.mean(temp2**2))
        
Uinf=np.mean(U[0,50,:])

#Calculating integral parameters
yy=np.flip(Y[:,loc,1])
uu=np.flip(Umean[:,loc])
y_corr = np.insert(yy,0,0,axis=0)
Vel = np.insert(uu,0,0,axis=0)
U_norm=np.zeros(len(Vel))
U_sub=np.zeros(len(Vel))
U_norm = Vel/Uinf 
U_sub=1-U_norm
disp_t=np.trapz(U_sub,y_corr)
mom_t=np.trapz(np.multiply(U_norm,U_sub),y_corr)
H=disp_t/mom_t




%matplotlib qt
plt.figure(1)
plt.plot(Ymean[:,loc]/disp_t,Urms[:,loc]/Uinf,'-o')
plt.plot(Ymean[:,loc]/disp_t,Urmsp[:,loc]/Uinf,'-d')
plt.plot(Ymean[:,loc]/disp_t,Urmsn[:,loc]/Uinf,'-s')
plt.grid()
plt.xlabel('Y/disp_t')
plt.ylabel('u_{rms}/Uinf')
plt.legend(['Total rms','+ve rms','-ve rms'])

plt.figure(2)
plt.rcParams['font.size'] = '16'
plt.subplot(5,1,1)
plt.scatter(Uflc[99,50,:]/Uinf,Vflc[99,50,:]/Uinf);plt.xlabel('u/Uinf');plt.ylabel('v/Uinf');plt.xlim([-0.1,0.1])
plt.subplot(5,1,2)
plt.scatter(Uflc[95,50,:]/Uinf,Vflc[95,50,:]/Uinf);plt.xlabel('u/Uinf');plt.ylabel('v/Uinf');plt.xlim([-0.1,0.1])
plt.subplot(5,1,3)
plt.scatter(Uflc[90,50,:]/Uinf,Vflc[90,50,:]/Uinf);plt.xlabel('u/Uinf');plt.ylabel('v/Uinf');plt.xlim([-0.1,0.1])
plt.subplot(5,1,4)
plt.scatter(Uflc[80,50,:]/Uinf,Vflc[80,50,:]/Uinf);plt.xlabel('u/Uinf');plt.ylabel('v/Uinf');plt.xlim([-0.1,0.1])
plt.subplot(5,1,5)
plt.scatter(Uflc[60,50,:]/Uinf,Vflc[60,50,:]/Uinf);plt.xlabel('u/Uinf');plt.ylabel('v/Uinf');plt.xlim([-0.1,0.1])




