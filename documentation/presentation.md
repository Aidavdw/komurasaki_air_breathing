# Introduction page
For calculating multi- cycle performance of microwave rocket, refilling performance is important.
Computational tools are important for designing to optimise this.
Changing positions of reed valves, etc.
Specific point of interest: new type: centre/cap valve, refilling in centre line
legacy code: old, limited, does not allow for this to be implemented

# Outline legacy code
## Physical:
- 2d, laminar flow model.
- Finite difference scheme, constant time stepping.
- 2017 for master thesis, focus on reed valves.
- Reed valves modeled as FEM beams, which based on their deflection add a source term.
- Reed valves use a lot of empirical parameters (obfuscated). Hard to generalise, easy to accidentally use wrongly.
## Meta
- input/output in (un-annotated) csv
- Changing configuration requires changing source C code, and recompiling. Goes against 'design' purpose, inhibiting rapid iteration.
- Extension is hard, old standard, poorly documented C code.
- A lot has been hard-coded.
- Many small mistakes: equation errors, memory leaks.
- Requires outdated compiling instruction set
- Basic parallel processing

# Outline new code
Long story short; had to extend, so instead I ended up re-building it from the start.
## Goals
- Easier to use in both analysis- and design problems: -> IO, documentation
- More freedom in changing of design: Allow valves on any wall, changing valve geometry, arbitrary domain shapes, other types of valves
- Easier to extend: Less hard coding, more generalised model -> modern design patterns, abstraction

# Introducing the new paradigm
- Complete Rewrite in modern C++ (2011+)
- Object-oriented design, placing design parameters into logical places and making extension easier for the large code base.
- Simplified build system for multiple platforms.
- Interface using a python library: Setting configuration using simple commands, getting output as numpy arrays and dataframes for quicker analysis.
- Detailed logging and error handling

# Fluid mechanics
2D compressible model in radial grid. Governing equations are relatively straightforward.
State variables: density, velocity (vector), internal energy.
## Viscosity, turbulence model
Investigation into effect of viscosity, and if turbulence models are required.
Simple comparative case has been set up in ANSYS CFX to compare the magnitude of its effect.
At the time scales and dimensions posed, less than 2% divergence.
Only when scaled up to 10x the size does it start exhibiting more of an effect (5% divergence).
For these stages of development, disregarding viscous stresses is sufficient, and turbulence can be disregarded.
# Solver
Numerical implementation has been optimised for multi-thread compatibility.
Solver uses a 4th order Runge Kutta scheme compared to the original, where operation can be done partially independently.
This allows for use of coarser time discretisation, resulting in overall faster computability.
Also lowers the numerical dissipation, which was relatively high in the legacy model. Exact magnitude has not been verified yet.
Flux splitting is unchanged
no shock fronts: AUSM-DV
shock fronts: HANEL
Sampled values are piece-wise reconstructed twice differentiable partial functions from MUSCL interpolation. This allows direct monitoring of numerical performance, allowing direct feedback to user if scheme requires altering.

# Program flow
Rough outline of how the program flows.
Understanding this, it can be divided up into different modules
Three main phases:
- set-up
- initialisation
- simulation
- decomission
## Set-up
Done through the API, user should set global variables and define what makes this case unique.
Allow creating fluid domain with boundary conditions. 
Defining where valves are located, and their configuration.
## Initialisation
Domain State: flow state intensive properties; density, velocity, internal energy.
Allow setting numbers from existing data, ambient conditions, or chapman-jouget theory for validation.
Valves State: deflection.
Allow manual setting, or 'matching' to the environment.
## Simulation
The main time loop, continues until the total time is simulated.
Because Runge-kutta iteration, extra internal loop.
Inside runge kutta loop: For the domains, calculate the fluxes, solving the euler equations.
Concurrently, update valve state, and calculate fluxes from this
Then combine these fluxes, and determine the accumulation using the fluid model.
Feed this value back, until runge-kutta iteration is finished. Then move to the next time step.
Note that valve & domain are drawn parallel.
Both write into their own buffer.
These buffers are later added together.
This latter operation is lighter, so heavier operation of calculation can be done concurrently.
## Decomissioning
Exporting all the data, cleaning up memory.

# Logical layout of the program
Up to here relatively similar to original program, but this had all these parameters hard coded in.
Adding new domains & valve config is hard.
Move away from monolythic -> modular
This means everything is easier to swap out, but also to wrap your head around.
Therefore, I will first outline the modular structure or logical layout.
I will move bottom-up; going from most global to most fundamental.
This will then be used to show the program flow later.
Not focus on individual variables, but rather focus on Object Ownership, relations, and data representation

## TwoDimensionalArray
Memory-safe way of representing a two-dimensional array of numbers with additional ghost cells around it.
Continuous in memory with low level access, but also friendly accessors
Prevent accidental usage:
    Setters and getters separated
    Indexing bound checking
    Separation between accessing 'main' and ghost

Analysis functions, used for debugging and runtime-checking
Empty checking, square, triangularity, symmetry, etc.
Matrix operations like decomposition, transposition.

## Field Quantity
Represents a specific flow variable field quantity.
e.g. Density, pressure, enthalpy.
Thin wrapper around TwoDimensionalArray.
Main purpose: multiple buffers.
Functionality for accessing at specific buffer, moving/copying buffers.

## Domain
Represents a physical area where flow can propogate.
Specifically, the geometry of it.
Internal properties implemented as FieldQuantities.
Makes it also easy to extend to include other flow properties if desired later.
### Dimensional data & Mapping Computation to physical space
Instead, changes the abstractly defined FieldQuantities to make them actually represent fluid parameters and properties.
Dimensionalise the grid into a physical cells with variable grid spacing.
Several variable grid spacings are available so that more resolution can be given to areas of higher importance.
Also lower resolution in the far field to allow dissipation at a more reasonable computational cost.
Old code only had constant spacing, but this generally not desirable to adequately model flow near boundaries.
Also works the other way around;
functionality to go from a physical place to a cell. This is useful for the valves.
### Implementation of fluid models
FieldQuantities represent intensive variables. With dimensionalisation, they can now be extensive.
This enables modeling of the actual flow, so this has been implemented at this level.
This Euler implementation will be discussed in more detail later, in the program flow.
### Boundary Conditions
Several types for several uses.
Can be individually set for each boundary.
As inviscid flow is assumed, Slip, No slip.
Both fully conserve energy, no flux.
Work by populating the ghost cells to virtually effect the flow.
### Connection to other domains.
Special case boundary condition for connection to other boundaries.
This uses the ghost cells to copy values from other domain.
This allows all domains to run independently from each-other.
Ghost cells themselves are discarded afterwards.

## SimCase - Manager
Because interfacing with python, easiest way to expose is with a manager-type object.
Main object an end-user will be using to add domains, valves, set their relations, and changing overall parameters. 
This is what SimCase is.
Invokes the main loop, giving an anchoring points for sub-models, domains & valves.
Simulation length and time step are stored individually because of how important they are.
### SolverSettings:
Change how solver works
- Tuning of numerical parameters
- Order of solving schemes (Runge-kutta)
### Ambient Conditions
Represents the 'ambient' outside of domains.
Is used to initialise some domains, and calculating values like exergy.
Completely static.
### RuntimeParameters
Meta-model.
time steps between exporting
preferences for logging

## Valves
Remember, extending valves was main driver of rewrite.
Main function of valves is allowing air from different domain to enter.
### Interface: child classing
This is functionality shared between all valves, regardless of internal work.
Allows for implementing Interface Pattern;
shared set of functionalities that abstracts the inner workings, and only exposes the main functionality to the main solver
This allows easy extension and creation of new types of valves.
Only properties that need to be known:
- which domains it connects
- where it is located (boundary & position)
- Hooks for updating valve state, setting initial conditions, and calculating/caching fluxes.
### Accumulating fluxes
Current way to accumulate fluxes is similar to original program.
Adds a source term source location.
This can span multiple cells. Actual implementation is dependent on the type of valve.
Unfortunately not fully physical.
Great discontinuities, lots of numerical loss.
As will become clear later with reed valve: Needs many empirical parameters.
#### Dynamic Boundary Condition flux model
Partially implemented: Dynamic boundary conditions.
Using ghost cells directly, allow flux between the two domains with a flow resistance terms. *source
Reduces empirical parameters drastically by using the dimensionality of the problem; Scales automatically.
## Reed Valve
Specific implementation of reed valve.
Existing Reed valve model is extended.
Modeled as a FEM beam.
Forces are based on sampled pressures around it.
Solved using cholesky decomposition & Newmark solving.
Assumed that only the deflection orthogonal to boundary is relevant, axial is discarded.
If verificiation shows this should be implemented this is relatively simple, but requires slight re-writing of sampling code.
State variables are therefore degrees of freedom, deflection, velocity.
Time integrated to set new deflections, keeping track of velocity.
Damping is also applied based on the properties of the fluid around it. This is done by re-using the calibrated C damping values from Florian (2017)
### Source term calculation
Source term dependent on deflection of tip only.
Distributed normalised to the scaled flux-through area (cross-sectional area)
## Cap Valve
Aforementioned new type of valve.
Unlike reed valve, only assumed to have 1 degree of freedom: axial deflection.
Therefore far fewer state variables.
Deflection modeled as dampened axial spring.
Damping is directly related to a directional drag coefficient; easier to close than to open.
Existing literature and experimental data can be used to calculate flow resistance terms. *source
Common process ^.
Based on geometry, aerodynamic properties. For novel configurations, these can be calculated with CFD or experimentally.

# Python Interface
Different topic; Right now looking mostly at internal.
Final goal however is doing analysis with it.
Therefore, time spent on making it not only extendable, but also usable.
As mentioned in introduction, python interface instead of CSV and config
This allows easier setting up, but also easier analysis.
Also easier integration into larger pipelines, such as using in conjunction with shock generation simulation.
## Setting up a case
## Running a case
## Analysing

# What still needs to be done
- Verification & Validation
- Finshing the reed valve implementation
- Updating documentation
- Set up case with cap valve
- Add Record system for valves & abstract classes
- dynamic time-stepping
- write solver for non-constant mesh spacing

# Any questions?