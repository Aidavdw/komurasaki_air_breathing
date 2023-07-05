#Added animation code to Aida's example_python.py

#modules for animation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from collections import defaultdict
import matplotlib.cm as cm
import sys
import os
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] ="C:/ffmpeg/bin/ffmpeg.exe"

# the pyd module is in a subdir, add here for convenience
from sys import version_info, path
path.insert(0, './build/Debug')


from datetime import datetime
from sys import version_info

if (version_info.major < 3):
    if (version_info.minor < 11):
        if (version_info.micro < 4):
            raise Exception("The binary requires at least python version 3.11.4")


# Hello kinoshita! Good luck :)

#Name of case and saved animation
MODE="Pressure"
title=MODE+".mp4"

# I made it all into one module! ^-^
import komurasakiairbreathing as ka

# Some values that are re-used.
length_of_tube = 0.5
height_of_tube = 0.028

# The main object that defines the 'simulation case'.
# Simcase(total simulation time, dt)
simcase = ka.SimCase()
simcase.dt = 10.0E-6
simcase.simulation_duration = 10.0E-4
nstep=int(simcase.simulation_duration/simcase.dt) #todo: A: add simple callback method in python library to get simcase.n_time_step

solver_settings = ka.SolverSettings() # for now, use the default solver settings. See cpp file for default values.
chapman_jouget = ka.ChapmanJougetInitialConditionParameters()
chapman_jouget.beam_power = 2000E3 # Watt
chapman_jouget.energy_absorption_coefficient = 1



# Mesh Spacing objects describe how he grid is laid out over the domain in one specific axis. Right now, only MeshSpacingType.constant is defined, and the others will probably crash the program.
# ka.MeshSpacing(
#   ka.MeshSpacingType,
#   length of the axis. Should be equal to the size that will be defined in AddDoman,
#   total amount of cells that will be distributed
#   characteristic density left (leave 0 for constant)
#   characteristic density right (leave 0 for constant)
# )

x_mesh_spacing = ka.MeshSpacing()
x_mesh_spacing.type = ka.MeshSpacingType.constant
x_mesh_spacing.amount_of_elements = 250
y_mesh_spacing = ka.MeshSpacing()
y_mesh_spacing.type = ka.MeshSpacingType.constant
y_mesh_spacing.amount_of_elements = 14

# A 2d, rectangular domain where the flow will be simulated
tube = simcase.AddDomain(
    1,                                                      # ID
    "Tube",                                                 # Name, used for outputting
    ka.Position(0,0),                                       # Where the bottom-left corner of this domain is placed
    (length_of_tube, height_of_tube),                       # The x,y dimensions of this domain
    (x_mesh_spacing, y_mesh_spacing),                       # Mesh spacing parameters; how the grid is distributed
    ka.InitialisationMethod.from_chapman_jouget_solution,   # How this domain is initialised
)

# Setting the boundaries for the domain. Should be pretty self-explanatory.
tube.SetBoundaryType(ka.Face.left, ka.BoundaryCondition.slip)
tube.SetBoundaryType(ka.Face.top, ka.BoundaryCondition.slip)
tube.SetBoundaryType(ka.Face.bottom, ka.BoundaryCondition.slip)

ambient = simcase.AddDomain(
    2,
    "Ambient",
    ka.Position(length_of_tube,0),
    (length_of_tube, height_of_tube),  
    (x_mesh_spacing, y_mesh_spacing),
    ka.InitialisationMethod.ambient_conditions,
)

ambient.SetBoundaryType(ka.Face.right, ka.BoundaryCondition.slip)
ambient.SetBoundaryType(ka.Face.top, ka.BoundaryCondition.slip)
ambient.SetBoundaryType(ka.Face.bottom, ka.BoundaryCondition.slip)

# the tube and the ambient share a boundary. This connects them.
simcase.ConnectBoundariesByName("Tube", ka.Face.right, "Ambient", ka.Face.left)

# Sets what we're interested in. After running the simulation, we can access these records.
recordName1="Tube_"+MODE
recordName2="Ambient_"+MODE
if MODE=="Density":
    simcase.AddRecord(tube.density.current_time_step,recordName1)
    simcase.AddRecord(ambient.density.current_time_step, recordName2)
    CBrabel="Density []"
elif MODE=="Temperature":
    simcase.AddRecord(tube.temperature.current_time_step,recordName1)
    simcase.AddRecord(ambient.temperature.current_time_step, recordName2)   
    CBrabel="Temperature []"
elif MODE=="Pressure":
    simcase.AddRecord(tube.pressure.current_time_step,recordName1)
    simcase.AddRecord(ambient.pressure.current_time_step, recordName2)   
    CBrabel="Presuure [Pa]"



# 数値魔法
start_time = datetime.now()
ka.DoSimulation(simcase)
end_time = datetime.now()

print("Simulation done! time ran: ", end_time - start_time)



#Animation
print("Nstep=",nstep)

#Get the range of variable in order to set the range of colorbar
var_min=0
var_max=0
for i in range(nstep):
    a=simcase.two_dimensional_array_records[recordName1].AsNumpyArray(i)
    b=simcase.two_dimensional_array_records[recordName2].AsNumpyArray(i)
    if var_min == 0:
        var_min = min(a.min(),b.min())
    else:
        if var_min > min(a.min(),b.min()):
            var_min = min(a.min(),b.min())
        if var_max < max(a.max(),b.max()):
            var_max = max(a.max(),b.max())
    

def animate(nt):
    fig.clear()
    #Create directory to store images
    exportPath = r'IMAGES'
    if not os.path.isdir(exportPath):
        os.makedirs(exportPath)
    
    #Display field in background
    a=simcase.two_dimensional_array_records[recordName1].AsNumpyArray(nt)
    b=simcase.two_dimensional_array_records[recordName2].AsNumpyArray(nt)
    xarr1=np.linspace(0,length_of_tube,x_mesh_spacing.amount_of_elements)
    xarr2=np.linspace(length_of_tube,length_of_tube*2,x_mesh_spacing.amount_of_elements)
    yarr=np.linspace(0,height_of_tube,y_mesh_spacing.amount_of_elements)
    plt.pcolormesh(xarr1,yarr,a,cmap='jet',vmin=var_min,vmax=var_max)
    plt.pcolormesh(xarr2,yarr,b,cmap='jet',vmin=var_min,vmax=var_max)

    cb = plt.colorbar(fraction=0.05,aspect=40,orientation="horizontal")
    cb.ax.tick_params(labelsize=20) 
    cb.set_label(CBrabel,fontweight='bold',fontsize=25)
	# cb.ax.xaxis.set_ticks_position('top')

    # Plot options
    fig.set_tight_layout(True)
    # plt.xlim(xmin*1000,xmax*1000)
    # plt.ylim(ymin*1000,ymax*1000)
    plt.xticks(color='k',size=20,fontweight='bold')
    plt.yticks(color='k',size=20,fontweight='bold')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(MODE+str(round(nt*simcase.dt*10**3,2))+' ms.',fontsize=35,fontweight='bold')
    plt.xlabel('Axial Direction (m)',fontsize=20,fontweight='bold')
    plt.ylabel('Radial Direction (m)',fontsize=20,fontweight='bold')
    
    # Show rocket walls
    wallthickness=0.002
    rectangle_x = [0,0,length_of_tube,length_of_tube]
    rectangle_y = [0,-wallthickness,-wallthickness,0]
    triangle_x = [-0.02,0,0]
    triangle_y = [0,height_of_tube+wallthickness,0]
    plt.fill([x for x in rectangle_x],[x for x in rectangle_y],'0.75',linewidth=0)
    plt.fill([x for x in triangle_x],[x for x in triangle_y],'0.75',linewidth=0)

    plt.draw()

    return
    #End of animate function

# EXPORT ANIMATION
fig, axes = plt.subplots(num=None, figsize=(15, 8), dpi=200, facecolor='w', edgecolor='k')
fig.clear()
plt.axis('equal')

anim=animation.FuncAnimation(fig,animate,range(nstep))
FFMpegWriter=animation.writers['ffmpeg']
writer=FFMpegWriter(fps=5, bitrate=100)
anim.save(title,writer=writer)


pause = input()