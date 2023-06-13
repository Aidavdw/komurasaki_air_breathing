# File Outline
 - ``extern/`` contains external dependencies. Right now, this is only [Pybind11](https://github.com/pybind/pybind11).
 - ``redundant/`` contains old versions of Florian's code, kept for reference.
 - ``main.cpp`` is a stand-alone version of the simulation program. This is only useful for development purposes, as the test case declaration is much more elaborate inside of CPP. Can be useful if you are extending the program.
 - ``pythonmodule.cpp`` is the main build file for the python module. It is built with the [cmake](https://cmake.org/) build configuration.
 - ``buildpythonmodule.bat`` is a shorthand way to build the python module on a windows system. It expects requires visual studio, and it requires msbuild to be added to path
 - ``private/`` and  ``public/`` contains the header files and source files respectively.
 - ``private/debug/`` and ``private/debug/`` contains files exclusive to debugging purposes.
 - ``private/pythoninterface/`` and ``private/pythoninterface/`` contains files exclusive to interfacing with python, what could not be put in pythonmodule.cpp.
 

 # Program flow
 The program essentially works in **four** steps:
 1. setting up the case. This is done in python using the API, or in main.cpp for debug builds. The main purpose of this is to set all the relevant parameters in a SimCase, define the domains, and register valves.
 2. Initialising the simulation using the provided information. This means calculating initial conditions, etc.
 3. Running the simulation~
 4. Processing of results (this is up to you!)

A graphical version of this overview is shown below.
 ![A diagram of the flow of the program](program_flow.svg)

## Setting up the simulation case
The ``SimCase`` object works as the master container. It contains global information such as time step, and total simulation duration. It orchestrates the simulation by driving domains and valves.
From it, you can create new **Domains** and **Valves**.
The ``Domain`` represents a domain where fluid is present. In the set-up period, you give it a size, location, and method to initialise with.

# Data Architecture
The simcase contains a master reference to all the domains and valves.
The domains have their properties described by a collection of FieldQuantities.
These FieldQuantities again have multiple buffers (values it can have simultaneously), implemented as TwoDimensionalArrays.

Valves use an interface architecture, so that the program only needs to know about its core functionality and not its implementation. This means that making a new type of valve is as simple as inheriting from *IValve*. Right now, only reed valves are implemented. They have their defomation implemented separately as FemDeformation.
![A diagram showing the data layout of the program](logical-layout.svg)
# Coding style

## Accessor consts
the TwoDimensionalArray, domains, and reed valves all have both const accessors and non-const accessors. read-only accessors are something like ``GetAt()``, whereas accessors are like ``At()``. This is done to preserve const-correctness, which makes sure you don't accidentally set a value instead of reading it. Try to maintain this if you are expanding the code, as it makes it easier to ensure you don't make these mistakes which are really hard to find!

# Build instructions
This program has 2 configurations;
- **standalone** is useful for extending the program as it keeps everything inside of c++. This means that it's much harder to read though. As I developed this on eclipse and visual studio, this configuration is accessed using ``air_breathing.sln``. Its main file is ``main.cpp``. This configuration has no external dependencies.
- **Python module** allows you to use the code through python. It is powered by [Pybind11](https://github.com/pybind/pybind11) and builds through [cmake](https://cmake.org/). its main file is ``pythonmodule.cpp``.

## Installing/updating Pybind11
Pybind11 is in this repository as a submodule. if it is missing, you must``git submodule init`` & ``git submodule update``.
## Python development headers
Compiling the python library requires you to have python development headers installed. If you don't, it might give you an error like *python.h not found*. These can be installed using *pip*, by ``pip install python-dev``
## Build flags
- ``_DEBUG`` adds a lot of run-time checking and outputs when it moves to the next operation. I still need to split these two functionalities.
- ``_CREATE_DUMP_FILES`` dumps all of the current_time_step values for every timestep to `debug_output/`. This is unformatted csv. Mostly useful if the python library breaks.
    
## Adding msbuild to path for use with buildpythonmodule.bat
using the ``buildpythonmodule.bat`` requires ``msbuild`` from visual studio to be added to path. In most normal installations, this program can be found at ``C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin``.