
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
XMAX = 1
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
        filename= case_dir + '/'+format(time[N_ITER_PLOT],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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


Pmin=p_min
Pmax=p_max
Umax=utot_max
Umin=utot_min
Tmin=t_min
Tmax=t_max
Rhomax=rho_max
Rhomin=rho_min

plt.close()

for l in range(0,N_image):



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
            filename= case_dir + '/'+format(time[N_ITER_PLOT+l],FORMAT)+'/'+str(k)+'/'+variable_array[nvar]+'.dat'
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




    # Plot of density, pressure and velocity fields
    plt.rc('text',usetex=True)
    plt.rc('font',family='times')




    if Quantity=='Pressure':
        for k in range(NSTART,NEND):
            plt.pcolormesh(xn[k],yn[k],p[k],cmap='jet',vmin=60000,vmax=150000)

    if Quantity=='Density': 
        for k in range(NSTART,NEND):
            plt.pcolormesh(xn[k],yn[k],rho[k],cmap='jet',vmin=0.5,vmax=1.5)

    if Quantity=='Temperature': 
        for k in range(NSTART,NEND):
            plt.pcolormesh(xn[k],yn[k],T[k],cmap='jet',vmin=300,vmax=1000)


    if Quantity=='Velocity': 
        for k in range(NSTART,NEND):
            if MACH_ON==1: 
                plt.pcolormesh(xn[k],yn[k],utot[k],cmap='jet',vmin=0.3,vmax=1.5)
            else: 
                plt.pcolormesh(xn[k],yn[k],utot[k],cmap='jet',vmin=0,vmax=1*340)

    # Axes of PLOT
    plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
    plt.tick_params(labelsize=20)
    plt.xticks(color='k',size=15,fontweight='bold')
    plt.yticks(color='k',size=15,fontweight='bold')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('Axial Direction (mm)',fontsize=15,fontweight='bold')
    plt.ylabel('Radial Direction (mm)',fontsize=15,fontweight='bold')


    for k in plot_stream:
        # xline = np.array([xn[k][i][0] for i in range(0,nxcell[k])])
        # yline = np.array([yn[k][0][i] for i in range(0,nycell[k])])
        # print(xc[k][:,:].T.shape)
        # print(u[k].shape)
        #plt.streamplot(xc[k][:,:].T,yc[k][:,:].T,u[k].T,v[k].T,color='k',density=(10,10),arrowsize=0.1,linewidth=0.1,arrowstyle='->')
         plt.quiver(xc[k][:,:].T,yc[k][:,:].T,u[k].T,v[k].T,color='k',scale=1,scale_units='inches',headwidth=0.05, linewidth=0.1)


    # Pressure contours (choose what figure and fields to plot)
    if plot_contours:
        plt.contour(xc[4][:,:].T,yc[4][:,:].T,p[4].T,300,colors='k',vmin=p_min,vmax=p_max+1)
        # plt.contour(xc[0][:,:].T,yc[0][:,:].T,p[0].T,300,colors='k',vmin=p_min,vmax=p_max+1)

    # Display reed valves
    if SHOW_VALVE==1:
        x_v = dict()
        y_v = dict()
        for k in range(n_valve):
            x_v[k] = np.zeros((n_fem+1,1),dtype=float)
            y_v[k] = np.zeros((n_fem+1,1),dtype=float)
            filename= case_dir + '/'+format(time[N_ITER_PLOT+l],FORMAT)+'/valve_'+ str(k) +'.dat'
            file = open(filename,"r")
            lineIndex=0;
            for line in file.readlines():
                data = [x for x in line.split(" ") if x!='']
                data = [x for x in line.split(" ") if x!='\n']
                data = [float(i) for i in data]
                x_v[k][lineIndex] = data[0];
                y_v[k][lineIndex] = data[1];
                lineIndex=lineIndex+1
            file.close()

    # Show mesh
    if PLOT_MESH==True:
        for n in range(4):
            plt.subplot(2,2,n+1)
            for k in range(Ndomain):
                plt.plot(xn[k],yn[k],color='k');
                plt.plot(xn[k].T,yn[k].T,color='k',linewidth=0.5);

    # Show plenum walls
    if SHOW_P_WALLS==True:
        xmin_1 = xn[4][-1][0]
        xmax_1 = xn[4][-1][0] + 0.002
        ymin_1 = yn[4][0][0]
        ymax_1 = yn[4][0][-1] + 0.001
        shape_x1 = [xmin_1,xmax_1,xmax_1,xmin_1,xmin_1]
        shape_y1 = [ymin_1,ymin_1,ymax_1,ymax_1,ymin_1]
        xmin_2 = xn[4][0][0]
        xmax_2 = xn[4][-1][0] + 0.002
        ymin_2 = yn[4][0][-1] -0.001
        ymax_2 = yn[4][0][-1] + 0.001
        shape_x2 = [xmin_2,xmax_2,xmax_2,xmin_2,xmin_2]
        shape_y2 = [ymin_2,ymin_2,ymax_2,ymax_2,ymin_2]
        
        plt.fill([x for x in shape_x1],[x for x in shape_y1],'0.75',linewidth=0)
        plt.fill([x for x in shape_x2],[x for x in shape_y2],'0.75',linewidth=0)

     # Show rocket walls
    if SHOW_R_WALLS==True:
        xmin_rec = 0
        xmax_rec = xn[0][-1][0]
        ymin_rec = yn[0][0][-1]
        ymax_rec = ymin_rec + 0.002
        rectangle_x = [xmin_rec,xmax_rec,xmax_rec,xmin_rec,xmin_rec]
        rectangle_y = [ymin_rec,ymin_rec,ymax_rec,ymax_rec,ymin_rec]
        xmin_tri = -0.1
        triangle_x = [xmin_tri,0,0,xmin_tri,xmin_tri]
        triangle_y = [0,0,ymax_rec,ymax_rec,0]
        
        plt.fill([x for x in rectangle_x],[x for x in rectangle_y],'0.75',linewidth=0)
        plt.fill([x for x in triangle_x],[x for x in triangle_y],'0.75',linewidth=0)


    # Show valves
    if SHOW_VALVE:
        for k in range(n_valve):
            plt.plot(x_v[k],y_v[k],linewidth=1,marker='o',color='w', markersize=1)

   
    if RESIZE == True :
        axes=plt.gca()
        axes.set_xlim([XMIN,XMAX])
        axes.set_ylim([YMIN,YMAX])

    if Quantity=='Pressure':
        plt.title('Pressure (Pa)',fontsize=15)
        plt.savefig('pressure'+str(N_ITER_PLOT+l), dpi=250)
    if Quantity=='Density':
        plt.title('Density ($kg/m^{3}$)',fontsize=15)
        plt.savefig('density'+str(N_ITER_PLOT+l), dpi=250)

    if Quantity=='Temperature':
        plt.title('Temperature ($K$)',fontsize=15)
        plt.savefig('temperature'+str(N_ITER_PLOT+l), dpi=250)

    if Quantity=='Velocity':
        if MACH_ON==1: 
            plt.title('Mach number',fontsize=15)
            plt.savefig('velocity'+str(N_ITER_PLOT+l), dpi=250)
        else: 
            plt.title('Velocity $||u||$ ($m/s$)',fontsize=15)
            plt.savefig('velocity'+str(N_ITER_PLOT+l), dpi=250)

    plt.close()




print('Start exporting video from generated field \n')
l=0
writergif = animation.PillowWriter(fps=30)
frames = [] # for storing the generated images
fig=plt.figure('movie')
for l in range(N_image): 
    if Quantity=='Pressure':
        img = plt.imread('pressure'+str(N_ITER_PLOT+l)+'.png')

    if Quantity=='Density': 
        img = plt.imread('density'+str(N_ITER_PLOT+l)+'.png')

    if Quantity=='Velocity': 
        img = plt.imread('velocity'+str(N_ITER_PLOT+l)+'.png')

    if Quantity=='Temperature': 
        img = plt.imread('temperature'+str(N_ITER_PLOT+l)+'.png')



    imgplot=plt.imshow(img, cmap='jet')

  
    frames.append([imgplot])

ani = animation.ArtistAnimation(fig, frames, interval=10, blit=True,repeat_delay=10)
 
if Quantity=='Pressure':
    ani.save('pressure_evolution.gif',writer=writergif)

if Quantity=='Density': 
    ani.save('density_evolution.gif',writer=writergif)

if Quantity=='Velocity': 
    ani.save('velocity_evolution.gif',writer=writergif)

if Quantity=='Temperature': 
    ani.save('temperature_evolution.gif',writer=writergif)


plt.show()
