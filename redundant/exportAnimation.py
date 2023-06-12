## Made by Florian (2017)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'


# USER PARAMETERS
PLOT_MESH = True           # To display the mesh grid above numerical results
EXPORT_PNG = True          # To export in .png format (as part of animation procedure)
EXPORT_EPS = False         # To export in .eps format (as part of animation procedure)
PLOT_STREAMLINES = False   # To plot streamlines on animation and images
PLOT_SOLIDS = True         # To display non-fluid regions on animations and images (and mesh display)
DISPLAY_MESH = True        # To display the mesh on a separate Python figure for checking
EXPORT_ANIMATION = True    # To export an animation
PLOT_VALVE = True          # Display reed valves
PLOT_P_WALLS = False        # Display plenum walls (only valid for the "plenum_rocket" case)
SHOW_R_WALLS = False        # Display rocket walls (only valid for the "plenum_rocket" case)


# NAME OF CASE AND OF SAVED ANIMATION
CASE='Laser'           # Name of case to load
MODE="rho"                       # Field to display
title=CASE+"_"+MODE+".mp4"        # Name of the final animation file


# CHOOSE BOUNDARIES OF THE DISPLAY WINDOW
USE_WINDOW = 1
XMIN = -0.1
YMIN = -0.01
XMAX = 0.55
YMAX = 0.15


# ROUNDING AXIS GRADING OF COLORBARS
rho_precision = 0.01
uv_precision = 0.1
p_precision = 1000
e_precision = 1000
T_precision =  1


# READING PARAMETERS FROM THE C# SIMULATION...
filename = CASE+'/Parameters.dat'
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
print(str(n_valve) + ' valve(s) of ' + str(n_fem) + ' elements each.')


# TIME VECTOR CORRESPONDING OF EXPORTED DATA
Tstart = 0
Nstep = int(1+Tsim/(printInterval*dt))
# Nstep = 30
time = np.zeros((Nstep),dtype=float)

for i in range(Nstep): 
    time[i] = Tstart+i*dt*printInterval


# MESH DICTIONARIES FOR MESH...
xc = dict()
yc = dict()
xn = dict()
yn = dict()

for x in range(Ndomain):
    xc[x] = 0
    yc[x] = 0
    xn[x] = 0
    yn[x] = 0


# ASSESSING SIZE OF COMPUTATIONAL DOMAIN
variable_array = ['xc']
Nvariable = len(variable_array)
Nx = np.zeros(Ndomain,dtype=int)
Ny = np.zeros(Ndomain,dtype=int)
for k in range(Ndomain):
    filename=CASE+"/"+format(time[0],FORMAT)+'/'+str(k)+'/'+variable_array[0]+'.dat'
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


# ALLOCATING MEMORY FOR MESH
for x in range(Ndomain):
    xn[x] = np.zeros((Nx[x]+1,Ny[x]+1),dtype=float)
    yn[x] = np.zeros((Nx[x]+1,Ny[x]+1),dtype=float)
    xc[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    yc[x] = np.zeros((Nx[x],Ny[x]),dtype=float)


# RECOVERING MESH
variable_array = ['x','y','xc','yc']
Nvariable = len(variable_array)
for k in range(Ndomain):
    for nvar in range(Nvariable):
        filename=CASE+"/"+format(time[0],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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


# PARAMETER DICTIONARIES
rho = dict()
u = dict()
v = dict()
p = dict()
e = dict()
T = dict()

for x in range(Ndomain):
    rho[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    u[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    v[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    p[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    e[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
    T[x] = np.zeros((Nx[x],Ny[x]),dtype=float)
variable_array = ['rho','u','v','p','E','T']
Nvariable = len(variable_array)


# LOOK FOR EXTREMUM FOR EACH FIELD VALUES OVER TIME
for nt in range(Nstep):
    for k in range(Ndomain):
        for nvar in range(Nvariable):
            filename=CASE+"/"+format(time[nt],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
            file = open(filename,"r")
            lineIndex = 0;
            for line in file.readlines():
                data = [x for x in line.split(" ") if x!='']
                data = [x for x in line.split(" ") if x!='\n']
                data = [float(i) for i in data]
                if nvar==0:
                    rho[k][lineIndex] = data
                    if nt==0 and k==0 and lineIndex==0:
                        rho_min = rho[k][0][0]
                        rho_max = rho_min
                    else:
                        rho_min = min(min(rho[k][lineIndex]),rho_min)
                        rho_max = max(max(rho[k][lineIndex]),rho_max)
                if nvar==1:
                    u[k][lineIndex] = data
                    if nt==0 and k==0 and lineIndex==0:
                        u_min = u[k][0][0]
                        u_max = u_min
                    else:
                        u_min = min(min(u[k][lineIndex]),u_min)
                        u_max = max(max(u[k][lineIndex]),u_max)
                if nvar==2:
                    v[k][lineIndex] = data
                    if nt==0 and k==0 and lineIndex==0:
                        v_min = v[k][0][0]
                        v_max = v_min
                    else:
                        u_min = min(min(v[k][lineIndex]),v_min)
                        u_max = max(max(v[k][lineIndex]),v_max)
                if nvar==3:
                    p[k][lineIndex] = data
                    if nt==0 and k==0 and lineIndex==0:
                        p_min = p[k][0][0]
                        p_max = p_min
                    else:
                        p_min = min(min(p[k][lineIndex]),p_min)
                        p_max = max(max(p[k][lineIndex]),p_max)
                if nvar==4:
                    e[k][lineIndex] = data
                    if nt==0 and k==0 and lineIndex==0:
                        e_min = e[k][0][0]
                        e_max = e_min
                    else:
                        e_min = min(min(e[k][lineIndex]),e_min)
                        e_max = max(max(e[k][lineIndex]),e_max)
                if nvar==5:
                    T[k][lineIndex] = data
                    if nt==0 and k==0 and lineIndex==0:
                        T_min = T[k][0][0]
                        T_max = T_min
                    else:
                        T_min = min(min(T[k][lineIndex]),T_min)
                        T_max = max(max(T[k][lineIndex]),T_max)
                lineIndex=lineIndex+1;
            file.close()


# SCALE EXTREMUM VALUES (SHOULD BE CUSTOMIZED BY USER)
rho_min = 0.2#round(rho_min/rho_precision)*rho_precision
rho_max = 1.4#math.ceil(rho_max/rho_precision)*rho_precision
u_min = round(u_min/uv_precision)*uv_precision
u_max = math.ceil(u_max/uv_precision)*uv_precision
v_min = round(v_min/uv_precision)*uv_precision
v_max = math.ceil(v_max/uv_precision)*uv_precision
p_min = round(p_min/p_precision)*p_precision
p_max = 300000##math.ceil(p_max/p_precision)*p_precision
e_min = round(e_min/e_precision)*e_precision
e_max = math.ceil(e_max/e_precision)*e_precision
T_min = round(T_min/T_precision)*T_precision
T_max = 800#math.ceil(T_max/T_precision)*T_precision


# CONFIRM EXTREMUM VALUES 
print("Min/Max density was found to be " + str(rho_min) + " / " + str(rho_max) + ".")
print("Min/Max U-velocity was found to be " + str(u_min) + " / " + str(u_max) + ".")
print("Min/Max V-velocity was found to be " + str(v_min) + " / " + str(v_max) + ".")
print("Min/Max pressure was found to be " + str(p_min) + " / " + str(p_max) + ".")
print("Min/Max energy was found to be " + str(e_min) + " / " + str(e_max) + ".")
print("Min/Max temperature was found to be " + str(T_min) + " / " + str(T_max) + ".")


# ALLOCATE MEMORY FOR VALVE DATA
xv = dict()
yv = dict()
for k in range(n_valve):
	xv[k] = np.zeros(n_fem+1,dtype=float)
	yv[k] = np.zeros(n_fem+1,dtype=float)


# ANIMATION FUNCTION (NEEDED TO CREATE IMAGES AND ANIMATIONS)...
def animate(nt):
    fig.clear();
    # Create directory to store images
    exportPath = r'IMAGES'
    if not os.path.isdir(exportPath):
        os.makedirs(exportPath)

    # Recover (rho,u,v,p,E) from output file
    for k in range(Ndomain):
        for nvar in range(Nvariable):
            filename=CASE+"/"+format(time[nt],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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
                if nvar==5:
                    T[k][lineIndex] = data
                lineIndex=lineIndex+1;
            file.close()

    # Recovering valve displacements
    for k in range(n_valve):
        filename=CASE+"/"+format(time[nt],FORMAT)+'/valve_'+str(k)+'.dat'
        file = open(filename,"r")
        lineIndex = 0;
        for line in file.readlines():
            data = [x for x in line.split(" ") if x!='']
            data = [x for x in line.split(" ") if x!='\n']
            data = [float(i) for i in data]
            xv[k][lineIndex] = data[0]
            yv[k][lineIndex] = data[1]
            lineIndex = lineIndex +1

    # Display field in background
    for k in range(Ndomain):
        if MODE=="rho":
            plt.pcolormesh(xn[k],yn[k],rho[k],cmap='jet',vmin=rho_min,vmax=rho_max)
        if MODE=="T" :
            plt.pcolormesh(xn[k],yn[k],T[k],cmap='jet',vmin=T_min,vmax=T_max)
        if MODE=="p" :
            plt.pcolormesh(xn[k],yn[k],p[k],cmap='jet',vmin=p_min,vmax=p_max)

    cb = plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
	# cb.ax.tick_params(labelsize=20) 
	# cb.set_label('Pressure (Pa)',fontweight='bold',fontsize=25)
	# cb.ax.xaxis.set_ticks_position('top')

    # Plot options
    fig.set_tight_layout(True)
    # plt.xlim(xmin*1000,xmax*1000)
    # plt.ylim(ymin*1000,ymax*1000)
    plt.xticks(color='k',size=20,fontweight='bold')
    plt.yticks(color='k',size=20,fontweight='bold')
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.title('Temperature field at time '+str(round(dt*printInterval*nt,2))+' s.',fontsize=35,fontweight='bold')
    plt.xlabel('Axial Direction (m)',fontsize=20,fontweight='bold')
    plt.ylabel('Radial Direction (m)',fontsize=20,fontweight='bold')
    if USE_WINDOW:
        axes=plt.gca()
        axes.set_xlim([XMIN,XMAX])
        axes.set_ylim([YMIN,YMAX])

    # Display valve distorsion
    if PLOT_VALVE:
        for k in range(n_valve):
    	   #plt.plot(xv[k].T,yv[k].T,linewidth=2,marker='o',color='w',markersize=0.01)
           plt.plot(xv[k],yv[k],linewidth=2,marker='o',color='w',markersize=3)


    # Show rocket walls
    if SHOW_R_WALLS==True:
        xmin_rec = 0
        xmax_rec = xn[0][-1][0]
        ymin_rec = yn[0][0][-1]
        ymax_rec = ymin_rec + 0.002
        rectangle_x = [xmin_rec,xmax_rec,xmax_rec,xmin_rec,xmin_rec]
        rectangle_y = [ymin_rec,ymin_rec,ymax_rec,ymax_rec,ymin_rec]
        xmin_tri = -0.05
        triangle_x = [xmin_tri,0,0,xmin_tri]
        triangle_y = [0,0,ymax_rec,0]
        plt.fill([x for x in rectangle_x],[x for x in rectangle_y],'0.75',linewidth=0)
        plt.fill([x for x in triangle_x],[x for x in triangle_y],'0.75',linewidth=0)


    # Show plenum walls
    if PLOT_P_WALLS==True:
        xmin_1 = xn[4][-1][0]
        xmax_1 = xn[4][-1][0] + 0.002
        ymin_1 = yn[4][0][0]
        ymax_1 = yn[4][0][-1] + 0.002
        shape_x1 = [xmin_1,xmax_1,xmax_1,xmin_1,xmin_1]
        shape_y1 = [ymin_1,ymin_1,ymax_1,ymax_1,ymin_1]
        xmin_2 = xn[4][0][0]
        xmax_2 = xn[4][-1][0] + 0.002
        ymin_2 = yn[4][0][-1] 
        ymax_2 = yn[4][0][-1] + 0.002
        shape_x2 = [xmin_2,xmax_2,xmax_2,xmin_2,xmin_2]
        shape_y2 = [ymin_2,ymin_2,ymax_2,ymax_2,ymin_2]
        plt.fill([x for x in shape_x1],[x for x in shape_y1],'0.75',linewidth=0)
        plt.fill([x for x in shape_x2],[x for x in shape_y2],'0.75',linewidth=0)


    # Update figure
    plt.draw()


    # Save figure into .png format (or .eps)
    if EXPORT_PNG:
        extension = '.png' 
        # exportFileName = exportPath+'/'+format(time[nt],FORMAT)+extension
        exportFileName = exportPath+'/'+str(nt)+extension
        print('Exporting :' + exportFileName)
        fig.savefig(exportFileName,dpi=500)
    if EXPORT_EPS:
        extension = '.eps' 
        exportFileName = exportPath+'/'+format(time[nt],FORMAT)+extension
        print('Exporting :' + exportFileName)
        fig.savefig(exportFileName,format='eps',dpi=500)

    return
    # END OF ANIMATE FUNCTION


# EXPORT ANIMATION
fig, axes = plt.subplots(num=None, figsize=(15, 8), dpi=50, facecolor='w', edgecolor='k')
fig.clear()
plt.axis('equal')
if EXPORT_ANIMATION:
    # CREATION OF ANIMATION FIGURE FOR DISPLAY
    anim = animation.FuncAnimation(fig, animate,Nstep,blit=True)
    FFMpegWriter = animation.writers['ffmpeg']
    writer = FFMpegWriter(fps=10, bitrate=100)
    anim.save(title,writer=writer)

# # END OF FILE
