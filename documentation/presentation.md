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

# Program flow
First, general steps the progam takes, not separated by components

# Logical layout of the program
Where variables are located
## SimCase - "Big Boss"

## Domain
Represents a physical area where flow can propogate.

## Field Quantity

## TwoDimensionalArray

## Valves
### Interface: child classing
### Reed Valve
### Cap Valve

# Program flow (logical)
The same program flow, but now separated by their logical partitioning in the code.

# Fluid mechanics

# Valve model
## Valve interface
not all of them- easy to extend.
## Reed valve re-implementation
## Cap Valve

# Python Interface
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