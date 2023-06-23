from datetime import datetime
from sys import version_info

if (version_info.major < 3):
    if (version_info.minor < 11):
        if (version_info.micro < 4):
            raise Exception("The binary requires at least python version 3.11.4")


# Hello kinoshita! Good luck :)

# I made it all into one module! ^-^
import komurasakiairbreathing as ka

# Some values that are re-used.
length_of_tube = 0.5
height_of_tube = 0.028

# The main object that defines the 'simulation case'.
# Simcase(total simulation time, dt)
simcase = ka.SimCase()
simcase.dt = 10.0E-7
simcase.simulation_duration = 10.0E-6

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
simcase.AddRecord(tube.density.current_time_step, "Tube_Density")
simcase.AddRecord(tube.pressure.current_time_step, "Tube_Pressure")




# 数値魔法
start_time = datetime.now()
ka.DoSimulation(simcase)
end_time = datetime.now()

print("Simulation done! time ran: ", end_time - start_time)

# Example of how to get the records at time step 5 for the TUbe Density record
a = simcase.two_dimensional_array_records['Tube_Density'].AsNumpyArray(5)
pause = input()