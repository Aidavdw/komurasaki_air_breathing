/*  This file contains functions that are related with Microwave Beaming and MSD theory:
    - eval_msd_function: Evaluates the value of function to solve to derive MSD Mach number.
    - eval_msd_deriv: Evaluates the value of the derivative of the former function.
    - solve_msd: Estimation of MSD Mach number based on Newton Raphson method. Returns the number of iterations required for sufficient convergence.

*/

#pragma once
#include "cell_values_container.h"

// forward declarations
class Domain;

// Represents a solution to a 1-d Chapman-Jouget detonation problem
struct ChapmanJougetDetonationSolution
{
    int iters_performed = 0;            // The amount of iterations performed solving the Chapman-Jouget equation before this solution was achieved.
    double m_msd = 0;                   // Mach number at which the detonation front propogates.
    double m1 = 0;                      // post-detonation mach number
    double detonation_velocity = 0;     // Velocity of the detonation front
    double l_exp = 0;                   // The position of the tail of the expansion region
    
    CellValues postExpansion;           // Post-expansion conditions; the 'plateau' region. This is what 'creeps in' from the rear, going towards the thrust wall.

    CellValues FieldPropertiesAtPosition(const double xPosition, const double gamma, const double specificGasConstant) const;
};

/* Function to solve to find MSD Mach number */
double eval_msd_function(double x, double c);

/* Derivative of function to solve to find MSD Mach number */
double eval_msd_deriv(double x, double c);

/* Compute initial pressure, density and temperature conditions based on Newton-Raphson's method. First return is the solution post-detonation, aka sol_1, and second return is the solution post-expansion, aka sol_2. */
ChapmanJougetDetonationSolution SolveChapmanJougetDetonationProblem(const double temperatureAmbient, const double pressureAmbient, const double eta, const double s0, const double idealGasConstant, const double gamma, const double lengthOfCombustionTube, const double radiusOfCombustionTube, const double convergenceThreshold = 1.0E-10);

// Populates a domain with the values in a Chapman-Jouget detonation. As the solution is 1D, the solution is the same across all the y-coordinates, and is only different for x coordinates.
void InitialiseDomainFromChapmanJougetDetonationSolution(Domain* domain, const ChapmanJougetDetonationSolution& sol, const double gamma, const double specificGasConstant);
