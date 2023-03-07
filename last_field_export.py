
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import math
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'




# USER PARAMETERS
N_ITER_PLOT = 0    # Iteration to plot (0 to 199, -1 to plot the last time step)
N_image=199              # Number of images wanted for the animated video
PLOT_MESH = False      # To display the mesh grid
SHOW_VALVE = True      # To show reed valves above solution fields
SHOW_R_WALLS = True   # To show rocket walls (valid for "plenum_rocket" case only)
SHOW_P_WALLS = False   # To show plenum walls (valid for "plenum rocket" case only)
NSTART = 0             # ID of first domain to plot (usually starts at 0)
NEND = 6            # ID of last index to plot (choose between 1 and the number of domains)
plot_stream = []       # ID of the domains for which to plot streamlines (ex: [0,1,2])
plot_contours = False  # To plot pressure contours above pressure field
stream_density = 1     # Density parameter when plotting streamlines
MACH_ON = False        # Plot total velocity in terms of Mach or absolute value


RESIZE=True
Quantity='Density'   #Temperature/Pressure/Velocity/Density 


# Boundaries of the display window of each subplot
USE_WINDOW = 1
XMIN = -0.1
YMIN = -0.01
XMAX = 0.2
YMAX = 0.1


# Name of case to display
case_dir = 'OUTPUT'


# Reading from parameter file of the designated case
filename = case_dir + '/Parameters.dat'
file = open(filename,"r")
Ndomain =int(file.readline())           # Number of domains
Tsim = float(file.readline())           # Total duration of computation (s)
dt = float(file.readline())             # Time step (m)
n_valve = int(file.readline())        # Number of reed valves
n_fem = int(file.readline())          # Number of finite elements per valve

line = file.readline()
xstart = [x for x in line.split(" ") if x!='']
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
nxcell = [int(x) for x in nxcell]

line = file.readline()
nycell = [x for x in line.split(" ") if x!='']
nycell = [x for x in line.split(" ") if x!='\n']
nycell = [int(x) for x in nycell]

dx = [xlength[i]/nxcell[i] for i in range(Ndomain)]
dy = [ylength[i]/nycell[i] for i in range(Ndomain)]

printInterval = int(file.readline())    # Number of time steps between two exports
FORMAT = str(file.readline()).strip("%")# Format for folder names
file.close()


# Display information to user
print('Number of domains: '+str(Ndomain))
print('Duration of simulation: '+str(Tsim))
print('Time step: '+str(dt))

print('XSTART: '+str(xstart))
print('YSTART: '+str(ystart))

print('XLENGTH: '+str(xlength))
print('YLENGTH: '+str(ylength))

print('Export Interval: '+str(printInterval)+' steps between exports ('+str(printInterval*dt*1000)+' ms)')
print('Writing format: '+FORMAT);
print('Current time of solution: '+str(N_ITER_PLOT*dt*printInterval))
print(str(n_valve) + ' valve(s) of ' + str(n_fem) + ' elements each.')





# Time vector of exported data
Tstart = 0
Nstep = int(1+Tsim/(float)(printInterval*dt))
time = np.zeros((Nstep),dtype=float)

for i in range(Nstep): 
    time[i] = Tstart+i*dt*printInterval

print('Total number of steps: ' + str(Nstep))


# Definition of dictionaries
xn = dict()
XN = dict()
yn = dict()
YN = dict()
xc = dict()
yc = dict()
XC = dict()
YC = dict()
rho = dict()
u = dict()
v = dict()
p = dict()
e = dict()
utot = dict()
T = dict()




# Import domain size information based on xc size
variable_array = ['xc']
Nvariable = len(variable_array)
Nx = np.zeros(Ndomain,dtype=int)
Ny = np.zeros(Ndomain ,dtype=int)
for k in range(Ndomain):
    filename=case_dir + '/'+format(time[N_ITER_PLOT],FORMAT)+'/'+str(k)+'/'+variable_array[0]+'.dat'
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

# Dictionaries to store solutions
for x in range(NSTART,NEND):
    xn[x] = np.zeros((Nx[x]+1,Ny[x]+1),dtype=float)
    yn[x] = np.zeros((Nx[x]+1,Ny[x]+1),dtype=float)
    xc[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    yc[x] = np.zeros((Nx[x],Ny[x]),dtype=float)

    rho[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    u[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    v[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    p[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    e[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    utot[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    T[x] = np.zeros((Nx[x],Ny[x]),dtype=float)

variable_array = ['x','y','xc','yc','rho','u','v','p','E','T']
Nvariable = len(variable_array)

  # Recover (rho,u,v,p,E) from output file
u_min=-1
u_max=1
v_min=-1
v_max=1
utot_max = 1
utot_min = 0
e_min=0
e_max=0




for k in range(NSTART,NEND):
    for nvar in range(Nvariable):
        filename= case_dir + '/'+format(time[N_ITER_PLOT+199],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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
            if nvar==4:
                rho[k][lineIndex] = data
                if lineIndex==0 and k==NSTART:
                    rho_min=min(data)
                    rho_max=max(data)
                else:
                    rho_min=min(min(data),rho_min)
                    rho_max=max(max(data),rho_max)
            if nvar==5:
                u[k][lineIndex] = data
                u_min=min(min(data),u_min)
                u_max=max(max(data),u_max)
            if nvar==6:
                v[k][lineIndex] = data
                v_min=min(min(data),v_min)
                v_max=max(max(data),v_max)
            if nvar==7:
                p[k][lineIndex] = data
                if lineIndex==0 and k==NSTART:
                    p_min=min(data)
                    p_max=max(data)
                else:
                    p_min=min(min(data),p_min)
                    p_max=max(max(data),p_max)
            if nvar==8:
                e[k][lineIndex] = data
                if lineIndex==0 and k==NSTART:
                    e_min=min(data)
                    e_max=max(data)
                else:
                    e_min=min(min(data),e_min)
                    e_max=max(max(data),e_max)
            if nvar==9:
                T[k][lineIndex] = data
                if lineIndex==0 and k==NSTART:
                    t_min=min(data)
                    t_max=max(data)
                else:
                    t_min=min(min(data),t_min)
                    t_max=max(max(data),t_max)
            lineIndex=lineIndex+1;
        file.close()


# Total velocity is deduced from U and V fields
for k in range(NSTART,NEND):
    for i in range(nxcell[k]):
        for j in range(nycell[k]):
            if MACH_ON==1:
                utot[k][i][j] = math.sqrt(u[k][i][j]*u[k][i][j] + v[k][i][j]*v[k][i][j])/math.sqrt(1.4*287*T[k][i][j])
                utot_max = max(utot_max,utot[k][i][j])
               
            else:
                utot[k][i][j] = math.sqrt(u[k][i][j]*u[k][i][j] + v[k][i][j]*v[k][i][j])
                utot_max = max(utot_max,utot[k][i][j])



np.savetxt('p17000_end.csv',p[0][:][:],delimiter=",")
np.savetxt('u17000_end.csv',u[0][:][:],delimiter=",")
np.savetxt('v17000_end.csv',v[0][:][:],delimiter=",")
np.savetxt('rho17000_end.csv',rho[0][:][:],delimiter=",")


