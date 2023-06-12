## TITLE: ANIMATED RESULTS OF THE C# MICROWAVE ROCKET COMPUTATION
## AUTHOR: FLORIAN NGUYEN, Jan. 2017
##
## Description: Reads raw .dat files in separate folders and
## regroups them for mesh plotting, creation of animations
## and images, etc.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import *
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'


# USER PARAMETERS
CASE = 'OUTPUT'
N_ITER_PLOT = 0     # Iteration to plot (-1 to plot the last time step)
plot_valve = True   # Display reed valves
plot_walls = True   # Display walls (not updated to take into account plenum and rocket)

# BOUNDARIES OF THE DISPLAYED WINDOW
USE_WINDOW = 0
XMIN = -0.1
YMIN = -0.001
XMAX = 0.55
YMAX = 0.06

# READING PARAMETERS FROM THE C# SIMULATION...
filename = CASE+'/Parameters.dat'
file = open(filename,"r")
Ndomain =int(file.readline())           # Number of domains
Tsim = float(file.readline())           # Total duration of computation (s)
dt = float(file.readline())             # Time step (m)
n_valve = int(file.readline())        # Number of reed valves
n_fem = int(file.readline())          # Number of finite elements per valve

line = file.readline()
xstart = [x for x in line.split(" ") if x!=b'']
xstart = [x for x in line.split(" ") if x!='\n']
xstart = [float(x) for x in xstart]

line = file.readline()
ystart = [x for x in line.split(" ") if x!='']
ystart = [x for x in line.split(" ") if x!='\n']
ystart = [float(x) for x in ystart]

line = file.readline()
xlength = [x for x in line.split(" ") if x!='']
xlength = [x for x in line.split(" ") if x!='\n']
xlength = [float(x) for x in xlength]

line = file.readline()
ylength = [x for x in line.split(" ") if x!='']
ylength = [x for x in line.split(" ") if x!='\n']
ylength = [float(x) for x in ylength]

line = file.readline()
nxcell = [x for x in line.split(" ") if x!='']
nxcell = [x for x in line.split(" ") if x!='\n']
nxcell = [float(x) for x in nxcell]

line = file.readline()
nycell = [x for x in line.split(" ") if x!='']
nycell = [x for x in line.split(" ") if x!='\n']
nycell = [float(x) for x in nycell]

printInterval = int(file.readline())    # Number of time steps between two exports
FORMAT = str(file.readline()).strip("%")# Format for folder names
file.close()


# DISPLAY INFORMATION...
print('Number of domains: '+str(Ndomain))
print('Duration of simulation: '+str(Tsim))
print('Time step: '+str(dt))

print('XSTART: '+str(xstart))
print('YSTART: '+str(ystart))

print('XLENGTH: '+str(xlength))
print('YLENGTH: '+str(ylength))

print('Export Interval: '+str(printInterval)+' steps between exports ('+str(printInterval*dt*1000)+' ms)')
print('Writing format: '+FORMAT);


# TIME VECTOR CORRESPONDING OF EXPORTED DATA
Tstart = 0
Nstep = int(1+Tsim/(printInterval*dt))
time = np.zeros((Nstep),dtype=float)

for i in range(Nstep): 
    time[i] = Tstart+i*dt*printInterval

# PARAMETER DICTIONARIES
xn = dict()
yn = dict()
xc = dict()
yc = dict()

# Reed valves
if plot_valve:
    x_v = dict()
    y_v = dict()
    for k in range(n_valve):
        x_v[k] = np.zeros((n_fem+1,1),dtype=float)
        y_v[k] = np.zeros((n_fem+1,1),dtype=float)
        filename=CASE+'/'+format(time[N_ITER_PLOT],FORMAT)+'/valve_'+ str(k) +'.dat'
        file = open(filename,"r")
        lineIndex=0;
        for line in file.readlines():
            data = [x for x in line.split(" ") if x!='']
            data = [x for x in line.split(" ") if x!='\n']
            data = [float(i) for i in data]
            # print(data)
            x_v[k][lineIndex] = data[0];
            y_v[k][lineIndex] = data[1];
            lineIndex=lineIndex+1
        file.close()

# Import domain size information based on xc size
variable_array = ['xc']
Nvariable = len(variable_array)
Nx = np.zeros(Ndomain,dtype=int)
Ny = np.zeros(Ndomain ,dtype=int)
for k in range(Ndomain):
    filename=CASE+'/'+format(time[N_ITER_PLOT],FORMAT)+'/'+str(k)+'/'+variable_array[0]+'.dat'
    file = open(filename,"r")
    lineIndex=0
    for line in file.readlines():
        data = [x for x in line.split(" ") if x!='']
        data = [x for x in line.split(" ") if x!='\n']
        data = [float(i) for i in data]
        Ny[k] = len(data)
        lineIndex=lineIndex+1
    Nx[k] = lineIndex
    file.close()

for x in range(Ndomain):
    xn[x] = np.zeros((Nx[x]+1,Ny[x]+1),dtype=float)
    yn[x] = np.zeros((Nx[x]+1,Ny[x]+1),dtype=float)
    xc[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    yc[x] = np.zeros((Nx[x],Ny[x]),dtype=float)

variable_array = ['x','y','xc','yc']
Nvariable = len(variable_array)

# Total dimensions
# NX = Nx[0]+Nx[1]+Nx[2]+Nx[3]
# NY = Ny[3]+Ny[4]
# XN = np.zeros((NX+1,NY+1),dtype=float)
# YN = np.zeros((NX+1,NY+1),dtype=float)
# XC = np.zeros((NX,NY),dtype=float)
# YC = np.zeros((NX,NY),dtype=float)

# Creation of a figure to display results (and PMAX to define max pressure over all time steps)
fig_mesh, axes = plt.subplots(num=None, figsize=(70, 8), dpi=50, facecolor='w', edgecolor='k')
fig_mesh.clear()
plt.axis('equal')

# Recover mesh from output file
for k in range(Ndomain):
    for nvar in range(Nvariable):
        filename=CASE+'/'+format(time[N_ITER_PLOT],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
        file = open(filename,"r")
        lineIndex = 0;
        for line in file.readlines():
            data = [x for x in line.split(" ") if x!='']
            data = [x for x in line.split(" ") if x!='\n']
            data = [float(i) for i in data]
            if nvar==0:
                xn[k][lineIndex] = data
            if nvar==1:
                yn[k][lineIndex] = data
            if nvar==2:
                xc[k][lineIndex] = data
            if nvar==3:
                yc[k][lineIndex] = data 
            lineIndex=lineIndex+1;
        file.close()

colors = cm.rainbow(np.linspace(0, 1, Ndomain))
for k in range(Ndomain):
    plt.plot(xn[k],yn[k],color=colors[k]);
    plt.plot(xn[k].T,yn[k].T,color=colors[k]);
    if plot_valve:
        for l in range(n_valve):
            plt.plot(x_v[l],y_v[l],color='k',marker='o',linewidth=0.5)

    plt.plot(xc[k],yc[k],color=colors[k],marker='+',linestyle='None');
    plt.plot(xc[k].T,yc[k].T,color=colors[k],marker='+',linestyle='None');

xmin_rec = 0
xmax_rec = xn[0][-1][0]
ymin_rec = yn[0][0][-1]
ymax_rec = ymin_rec + 0.002
rectangle_x = [xmin_rec,xmax_rec,xmax_rec,xmin_rec,xmin_rec]
rectangle_y = [ymin_rec,ymin_rec,ymax_rec,ymax_rec,ymin_rec]
plt.fill([x for x in rectangle_x],[x for x in rectangle_y],'0.75',linewidth=0)

xmin_tri = -0.05
triangle_x = [xmin_tri,0,0,xmin_tri]
triangle_y = [0,0,ymax_rec,0]
plt.fill([x for x in triangle_x],[x for x in triangle_y],'0.75',linewidth=0)

if USE_WINDOW==1:
    axes=plt.gca()
    axes.set_xlim([XMIN,XMAX])
    axes.set_ylim([YMIN,YMAX])

plt.rc('text',usetex=True)
plt.rc('font',family='times')
plt.xlabel('Axial Direction (mm)',fontsize=25,fontweight='bold')
plt.ylabel('Radial Direction (mm)',fontsize=25,fontweight='bold')
plt.xticks(color='k',size=25,fontweight='bold')
plt.yticks(color='k',size=25,fontweight='bold')

# Set limits on axes
# axes=plt.gca()
# axes.set_xlim([-0.01,0.45])
# axes.set_ylim([-0.001,0.061])

# if USE_WINDOW:
#     # axes.set_autoscale_on(False) 
#     # axes.set_xbound(lower=XMIN,upper=XMAX)
#     axes.set_ybound(lower=YMIN,upper=YMAX)

# Update figure
plt.show()

# END OF FILE