from datetime import datetime

# I made it all into one module! ^-^
import komurasakiairbreathing as ka

# Some values that are re-used.
length_of_tube = 0.5
height_of_tube = 0.028

# The main object that defines the 'simulation case'.
# Simcase(total simulation time, dt)
simcase = ka.SimCase(0.1, 10.0E-7)

# Mesh Spacing objects describe how he grid is laid out over the domain in one specific axis. Right now, only MeshSpacingType.constant is defined, and the others will probably crash the program.
# ka.MeshSpacing(
#   ka.MeshSpacingType,
#   length of the axis. Should be equal to the size that will be defined in AddDoman,
#   total amount of cells that will be distributed
#   characteristic density left (leave 0 for constant)
#   characteristic density right (leave 0 for constant)
# )
x_mesh_spacing = ka.MeshSpacing(ka.MeshSpacingType.constant, length_of_tube, 250, 0, 0)
y_mesh_spacing = ka.MeshSpacing(ka.MeshSpacingType.constant, height_of_tube, 14, 0, 0)

# A 2d, rectangular domain where the flow will be simulated
tube = simcase.AddDomain(
    1,                                                      # ID
    "Tube",                                                 # Name, used for outputting
    ka.Position(0,0),                                       # Where the bottom-left corner of this domain is placed
    (length_of_tube, height_of_tube),                       # The x,y dimensions of this domain
    (x_mesh_spacing, y_mesh_spacing),                       # Mesh spacing parameters; how the grid is distributed
    ka.InitialisationMethod.from_chapman_jouget_solution,   # How this domain is initialised
    2,                                                      # amount of ghost cells. Leave at 2 for now, will remove soon.
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
    2
)

ambient.SetBoundaryType(ka.Face.right, ka.BoundaryCondition.slip)
ambient.SetBoundaryType(ka.Face.top, ka.BoundaryCondition.slip)
ambient.SetBoundaryType(ka.Face.bottom, ka.BoundaryCondition.slip)

# the tube and the ambient share a boundary. This connects them.
simcase.ConnectBoundariesByName("Tube", ka.Face.right, "Ambient", ka.Face.left)

# Sets what we're interested in. After running the simulation, we can access these records.
simcase.AddRecord(tube.density.current_time_step, "Tube_Density")



# 数値魔法
start_time = datetime.now()
ka.DoSimulation(simcase)
end_time = datetime.now()

print("Simulation done! time ran: ", end_time - start_time)
pause = input()