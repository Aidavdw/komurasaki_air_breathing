## TITLE: EXTRACT DATA FROM A SPECIFIED LINE AND EXPORT IT
## FOR LATER COMPARISON WITH OPENFOAM DATA ON EXCEL
## AUTHOR: FLORIAN NGUYEN, Jan. 2017
##
## Description: Read data from the last time step, identify
## the indexes corresponding to the area of interest (line
## to be specified by the user) and extract in a .txt file
## that can later be used for comparison with OpenFOAM data
## on Excel.


##
## HEADER
##

import numpy as np
import matplotlib.pyplot as plt
from math import *
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation

##
## END OF HEADER
##


# USER PARAMETERS
N_ITER_PLOT = -1          # Iteration to plot (-1 to plot the last time step)
x_line = 1
y_line_start = 0.2
y_line_end = 1

# OPENFOAM PARAMETERS (needed to correct rho when exporting)
R_gas = 8.3144598
NORM_MOLAR_WEIGHT = 11.64 # (kg/mol)
CODE_MOLAR_WEIGHT = 0.02896 # (kg/mol)
RHO_RATIO = NORM_MOLAR_WEIGHT/CODE_MOLAR_WEIGHT


# READING PARAMETERS FROM THE C# SIMULATION...
filename = 'OUTPUT/Parameters.dat'
file = open('OUTPUT/Parameters.dat',"r")
Ndomain =int(file.readline()) # Number of domains
Tsim = float(file.readline())
dt = float(file.readline())
dx = float(file.readline())
dy = float(file.readline())
printInterval = int(file.readline())
FORMAT = str(file.readline()).strip("%")
file.close()

# Tsim = 3

# DISPLAY INFORMATION...
print('Number of domains: '+str(Ndomain))
print('Duration of simulation: '+str(Tsim))
print('Time step: '+str(dt))
print('DX: '+str(dx*1000)+' mm ; DY: '+str(dy*1000)+' mm')
print('Export Interval: '+str(printInterval)+' steps between exports ('+str(printInterval*dt*1000)+' ms)')
print('Writing format: '+FORMAT);


# TIME VECTOR CORRESPONDING OF EXPORTED DATA
Tstart = 0
Nstep = int(1+Tsim/(printInterval*dt))
time = np.zeros((Nstep),dtype=float)
for i in range(Nstep): 
    time[i] = Tstart+i*dt*printInterval


# MESH DICTIONARIES FOR MESH...
xc = dict()
yc = dict()
XC = dict()
YC = dict()

for x in range(Ndomain):
    xc[x] = 0
    yc[x] = 0
    XC[x] = 0
    YC[x] = 0


# GRID INFORMATION
variable_array = ['xc', 'yc']
Nvariable = len(variable_array)
Nx = np.zeros(Ndomain,dtype=int)
Ny = np.zeros(Ndomain ,dtype=int)
for k in range(Ndomain):
    for nvar in range(Nvariable):
        filename='OUTPUT/4.850000/'+str(k)+'/'+variable_array[nvar]+'.dat'
        filename='OUTPUT/'+format(time[N_ITER_PLOT],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
        file = open(filename,"r")
        for line in file.readlines():
            data = [x for x in line.split(" ") if x!='']
            data = [float(i) for i in data]
        if nvar==0:
            Nx[k] = len(data)
            xc[k] = data
        if nvar==1:
            Ny[k] = len(data)
            yc[k] = data


# PARAMETER DICTIONARIES
rho = dict()
u = dict()
v = dict()
p = dict()
e = dict()

for x in range(Ndomain):
    rho[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    u[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    v[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    p[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    e[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
variable_array = ['rho','u','v','p','E']
Nvariable = len(variable_array)

# Recover (rho,u,v,p,E) from output file
for k in range(Ndomain):
    for nvar in range(Nvariable):
        # filename='OUTPUT/4.850000/'+str(k)+'/'+variable_array[nvar]+'.dat'
        filename='OUTPUT/'+format(time[N_ITER_PLOT],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
        file = open(filename,"r")
        lineIndex = 0;
        for line in file.readlines():
            data = [x for x in line.split(" ") if x!='']
            data = [x for x in line.split(" ") if x!='\n']
            data = [float(i) for i in data]
            if nvar==0:
                rho[k][lineIndex] = data
            if nvar==1:
                u[k][lineIndex] = data
            if nvar==2:
                v[k][lineIndex] = data
            if nvar==3:
                p[k][lineIndex] = data
            if nvar==4:
                e[k][lineIndex] = data
            lineIndex=lineIndex+1;

# Merge all domains into one
NX = Nx[0]+Nx[2]
NY = Ny[0]+Ny[2]
xstart = [0,0,Nx[0]]
ystart = [0,Ny[0],Ny[0]]
X = xc[0]+xc[2]
Y = yc[0]+yc[2]
RHO = np.zeros((NX,NY),dtype=float)
U = np.zeros((NX,NY),dtype=float)
V = np.zeros((NX,NY),dtype=float)
P = np.zeros((NX,NY),dtype=float)
E = np.zeros((NX,NY),dtype=float)
VTOT = np.zeros((NX,NY),dtype=float)
MTOT = np.zeros((NX,NY),dtype=float)
T = np.zeros((NX,NY),dtype=float)
for x in range(Ndomain):
    for i in range(Nx[x]):
        for j in range(Ny[x]):
            RHO[i+xstart[x]][j+ystart[x]] = rho[x][i][j]
            U[i+xstart[x]][j+ystart[x]] = u[x][i][j]
            V[i+xstart[x]][j+ystart[x]] = v[x][i][j]
            P[i+xstart[x]][j+ystart[x]] = p[x][i][j]
            E[i+xstart[x]][j+ystart[x]] = e[x][i][j]
            T[i+xstart[x]][j+ystart[x]] = p[x][i][j]/(rho[x][i][j]*287)
            VTOT[i+xstart[x]][j+ystart[x]] = np.sqrt(u[x][i][j]*u[x][i][j]+v[x][i][j]*v[x][i][j])
            MTOT[i+xstart[x]][j+ystart[x]] = np.sqrt(u[x][i][j]*u[x][i][j]+v[x][i][j]*v[x][i][j])/(sqrt(1.4*p[x][i][j]/rho[x][i][j]))

def find_nearest(array,value):
    index=0
    nearVal = array[0]
    for i in range(len(array)):
        if abs(array[i]-value)<=abs(nearVal-value):
            index = i
            nearVal = array[i]
    return index

## Locate index of interest
nx_line_prev = find_nearest(X,x_line)
nx_line_next = nx_line_prev + 1
ny_line_start = find_nearest(Y,y_line_start)
ny_line_end = find_nearest(Y,y_line_end)

print("x = " + str(x_line) + "m is fitted at " + str(X[nx_line_prev]) + "m (next at " + str(X[nx_line_next]) + "m).")
print("ymin = " + str(y_line_start) + "m is fitted at " + str(Y[ny_line_start]) + "m and ymax = " + str(y_line_end) + "m at " + str(Y[ny_line_end]) + "m.")

## Open file with given filename: export mesh and parameters
filename='48x16_MESH_PREV.dat'
file = open(filename,"w")

## Write data
a = 0
# for j in range(0,ny_line_end - ny_line_start + 1):
for i in range(ny_line_start,ny_line_end + 1): 
    # i = ny_line_end - j
    # a = sqrt(1.4*287*T[nx_line_prev][i])
    # file.write(str(X[nx_line_prev]) + " " + str(Y[i]) + " " + str(RHO[nx_line_prev][i]*RHO_RATIO) + " " + str(P[nx_line_prev][i]) + " " + str(T[nx_line_prev][i]) + " " + str(U[nx_line_prev][i]/a) + " " + str(V[nx_line_prev][i])+ " ")
    # a = sqrt(1.4*287*T[nx_line_next][i])
    # file.write(str(X[nx_line_next]) + " " + str(Y[i]) + " " + str(RHO[nx_line_next][i]*RHO_RATIO) + " " + str(P[nx_line_next][i]) + " " + str(T[nx_line_next][i]) + " " + str(U[nx_line_next][i]/a) + " " + str(V[nx_line_next][i]) + "\n")
    file.write(str(X[nx_line_prev]) + " " + str(Y[i]) + " " + str(RHO[nx_line_prev][i]) + " " + str(P[nx_line_prev][i]) + " " + str(T[nx_line_prev][i]) + " " + str(U[nx_line_prev][i]) + " " + str(V[nx_line_prev][i]) + " ")
    file.write(str(X[nx_line_next]) + " " + str(Y[i]) + " " + str(RHO[nx_line_next][i]) + " " + str(P[nx_line_next][i]) + " " + str(T[nx_line_next][i]) + " " + str(U[nx_line_next][i]) + " " + str(V[nx_line_next][i]) + str(U[NX][i]*rho[NX][i]) + "\n")
file.close()


# # END OF FILE