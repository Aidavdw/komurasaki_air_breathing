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
    ChapmanJougetDetonationSolution(const CellValues& postExpansion) :
        postExpansion(postExpansion)
    {

    }

    int iters_performed = 0;            // The amount of iterations performed solving the Chapman-Jouget equation before this solution was achieved.

    double m_msd = 0;                   // Mach number at which the detonation front propogates.
    double m1 = 0;                      // post-detonation mach number
    double detonation_velocity = 0;     // Velocity of the detonation front

    //CellValues postDetonation;          // Post-detonation conditions, immediately behind the detonation front. This is what holds at the thrust wall. (subscript 1). Not serialised as of now.
    CellValues postExpansion;           // Post-expansion conditions. This is what 'creeps in' from the rear, going towards the thrust wall.

    //double p1 = 0;                      // post-detonation pressure
    //double u1 = 0;                      // post-detonation speed ? not sure what it represents
    //double rho1 = 0;                    // post-detonation density

    //double p2 = 0;                      // post-expansion pressure
    //double rho2 = 0;                    // post-expansion density
    double l_exp = 0;                   // The position of the tail of the expansion region

    CellValues FieldPropertiesAtPosition(const double xPosition, const double gamma, const double tubeRadius) const;
};

/* Function to solve to find MSD Mach number */
double eval_msd_function(double x, double c);

/* Derivative of function to solve to find MSD Mach number */
double eval_msd_deriv(double x, double c);

//todo: rename R0 parameter, as its confusing.
/* Compute initial pressure, density and temperature conditions based on Newton-Raphson's method.  */
ChapmanJougetDetonationSolution SolveChapmanJougetDetonationProblem(const double temperatureAmbient, const double pressureAmbient, const double ETA, const double S0, const double R, const double GAMMA, const double lengthOfCombustionTube, const double RadiusOfCombustionTube, const double convergenceThreshold = 1.0E-10);

// Populates a domain with the values in a Chapman-Jouget detonation. As the solution is 1D, the solution is the same across all the y-coordinates, and is only different for x coordinates.
void InitialiseDomainFromChapmanJougetDetonationSolution(Domain* domain, const ChapmanJougetDetonationSolution& sol, const double gamma);
