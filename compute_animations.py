## TITLE: ANIMATED RESULTS OF THE C# MICROWAVE ROCKET COMPUTATION
## AUTHOR: FLORIAN NGUYEN, Jan. 2017
##
## Description: Reads raw .dat files in separate folders and
## regroups them for mesh plotting, creation of animations
## and images, etc.


##
## HEADER
##

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import floor
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation

##
## END OF HEADER
##


# USER PARAMETERS
PLOT_MESH = True           # To display the mesh grid
EXPORT_PNG = True          # To export in .png format (as part of animation procedure)
EXPORT_EPS = False         # To export in .eps format (as part of animation procedure)
PLOT_STREAMLINES = False   # To plot streamlines on animation and images
PLOT_SOLIDS = True         # To display non-fluid regions on animations and images (and mesh display)
DISPLAY_MESH = True        # To display the mesh as Python figure
EXPORT_ANIMATION = True    # To export an animation

# READING PARAMETERS FROM THE C# SIMULATION...
filename = 'OUTPUT/Parameters.dat'

file = open(filename,"r")

Ndomain =int(file.readline())           # Number of domains
Tsim = float(file.readline())           # Total duration of computation (s)
dt = float(file.readline())             # Time step (m)
dx = float(file.readline())             # Grid size (X)
dy = float(file.readline())             # Grid size (Y)
printInterval = int(file.readline())  # Number of time steps between two exports
FORMAT = str(file.readline()).strip("%")# Format for folder names
file.close()


# To limit the number of time steps to include, manually modify Tsim
# Tsim = 3.5


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

print("Time steps to be exported: "+ str(time))

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

# Import grid information
variable_array = ['xc', 'yc']
Nvariable = len(variable_array)
Nx = np.zeros(Ndomain,dtype=int)
Ny = np.zeros(Ndomain ,dtype=int)
for nt in range(Nstep):
    for k in range(Ndomain):
        for nvar in range(Nvariable):
            filename='OUTPUT/'+format(time[nt],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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
            file.close()


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

# Creation of a figure to display results (and PMAX to define max pressure over all time steps)
fig, axes = plt.subplots(num=None, figsize=(15, 10), dpi=50, facecolor='w', edgecolor='k')
PMAX = 0

# Look for the maximal pressure recorded during all time steps
for nt in range(Nstep):
    for k in range(Ndomain):
        for nvar in range(Nvariable):
            filename='OUTPUT/'+format(time[nt],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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
            file.close()

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

# Confirm value of max pressure
PMAX = floor(PMAX)
if (PMAX % 10):
    PMAX = PMAX + (10 - PMAX % 10)
print("Max pressure was found to be " + str(PMAX))

# Recover (rho,u,v,p,E) from output file
for k in range(Ndomain):
    for nvar in range(Nvariable):
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
# Xmin/max and Ymin/max of the figure
xmin = xc[0][0]-dx/2

# xmax = 2*xc[0][Nx[0]-1] # To partly remove outlet region
xmax = xc[2][Nx[2]-1]-dx/2
ymin = yc[0][0]-dy/2
ymax = yc[2][Ny[2]-1]-dy/2

# Coordinates are converted so that data is plotted at cell centres
[Xg,Yg] = np.meshgrid(X,Y)
Xg=Xg.T-dx/2
Yg=Yg.T-dy/2

# Obstacle definition
if PLOT_SOLIDS:
    xmin_s = xc[2][0]-dx/2
    xmax_s = xc[2][Nx[2]-1]+dx/2
    ymin_s = yc[0][0]-dy/2
    
    ymax_s = yc[0][Ny[0]-1]+dy/2
    rectangle_x = [xmin_s,xmax_s,xmax_s,xmin_s,xmin_s]
    rectangle_y = [ymin_s,ymin_s,ymax_s,ymax_s,ymin_s]

# CREATION OF FIGURE TO DISPLAY RESULTS
fig, axes = plt.subplots(2,1, figsize=(34, 18), dpi=50, facecolor='w', edgecolor='k')

# Subplot 1: Pressure
plt.subplot(2,2,1)
plt.pcolormesh(Xg*1000,Yg*1000,P,cmap='jet')
cb1 = plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
# cb1.set_label('Pressure (Pa)',fontweight='bold',fontsize=20)
# cb1.ax.xaxis.set_ticks_position('top')
cb1.ax.tick_params(labelsize=20) 
plt.xticks(color='k',size=15,fontweight='bold')
plt.yticks(color='k',size=15,fontweight='bold')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Axial Direction (mm)',fontsize=20,fontweight='bold')
plt.ylabel('Radial Direction (mm)',fontsize=20,fontweight='bold')
plt.title('Pressure (Pa)',fontsize=30,fontweight='bold')
plt.xlim(xmin*1000,xmax*1000)
plt.ylim(ymin*1000,ymax*1000)
plt.fill([1000*x for x in rectangle_x],[1000*x for x in rectangle_y],'0.75')

#Subplot 2: Velocity OR Mach number
plt.subplot(2,2,2)
plt.pcolormesh(Xg*1000,Yg*1000,MTOT,cmap='jet')
cb2 = plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
# cb2.set_label('Mach number',fontweight='bold',fontsize=20)
# cb2.ax.xaxis.set_ticks_position('top')
cb2.ax.tick_params(labelsize=20) 
plt.xticks(color='k',size=15,fontweight='bold')
plt.yticks(color='k',size=15,fontweight='bold')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Axial Direction (mm)',fontsize=20,fontweight='bold')
plt.ylabel('Radial Direction (mm)',fontsize=20,fontweight='bold')
plt.title('Mach number',fontsize=30,fontweight='bold')
plt.xlim(xmin*1000,xmax*1000)
plt.ylim(ymin*1000,ymax*1000)
plt.fill([1000*x for x in rectangle_x],[1000*x for x in rectangle_y],'0.75')

# Subplot 3: Streamlines and contours
plt.subplot(2,2,3)
plt.contourf(Xg*1000,Yg*1000,VTOT)
cb3 = plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
# cb3.set_label('Absolute speed (m/s)',fontweight='bold',fontsize=20)
# cb3.ax.xaxis.set_ticks_position('top')
cb3.ax.tick_params(labelsize=20) 
plt.streamplot(Xg.T*1000,Yg.T*1000,U.T,V.T,color='k',density=2,arrowsize=1,linewidth=0.8,arrowstyle='->')
plt.xticks(color='k',size=15,fontweight='bold')
plt.yticks(color='k',size=15,fontweight='bold')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Axial Direction (mm)',fontsize=20,fontweight='bold')
plt.ylabel('Radial Direction (mm)',fontsize=20,fontweight='bold')
plt.title('Absolute velocity (m/s)',fontsize=30,fontweight='bold')
plt.xlim(xmin*1000,xmax*1000)
plt.ylim(ymin*1000,ymax*1000)
plt.fill([1000*x for x in rectangle_x],[1000*x for x in rectangle_y],'0.75')

# Subplot 4: Temperature
plt.subplot(2,2,4)
plt.pcolormesh(Xg*1000,Yg*1000,T,cmap='jet')
cb1 = plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
# cb1.set_label('Pressure (Pa)',fontweight='bold',fontsize=20)
# cb1.ax.xaxis.set_ticks_position('top')
cb1.ax.tick_params(labelsize=20) 
plt.xticks(color='k',size=15,fontweight='bold')
plt.yticks(color='k',size=15,fontweight='bold')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Axial Direction (mm)',fontsize=20,fontweight='bold')
plt.ylabel('Radial Direction (mm)',fontsize=20,fontweight='bold')
plt.title('Temperature (K)',fontsize=30,fontweight='bold')
plt.xlim(xmin*1000,xmax*1000)
plt.ylim(ymin*1000,ymax*1000)
plt.fill([1000*x for x in rectangle_x],[1000*x for x in rectangle_y],'0.75')

# Title of the figure
fig.suptitle('Simulation results at time '+str('{0:.2f}'.format(dt*printInterval*time[N_ITER_PLOT]*1000))+' ms.',fontsize=35,fontweight='bold')
fig.set_tight_layout(True)
plt.show()
