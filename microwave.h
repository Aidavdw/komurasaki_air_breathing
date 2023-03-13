#include "parameters.h"

/*  This file contains functions that are related with Microwave Beaming and MSD theory:
    - eval_msd_function: Evaluates the value of function to solve to derive MSD Mach number.
    - eval_msd_deriv: Evaluates the value of the derivative of the former function.
    - solve_msd: Estimation of MSD Mach number based on Newton Raphson method. Returns the number of iterations required for sufficient convergence.

*/

// Represents a solution to a 1-d Chapman-Jouget detonation problem
struct ChapmanJougetDetonationSolution
{
    ChapmanJougetDetonationSolution();

    int iters_performed = 0;

    double m_msd = 0;
    double p1 = 0;
    double u1 = 0;
    double rho1 = 0;
    double m1 = 0;
    double p2 = 0;
    double rho2 = 0;
    double l_exp = 0;
};

/* Function to solve to find MSD Mach number */
double eval_msd_function(double x, double c);

/* Derivative of function to solve to find MSD Mach number */
double eval_msd_deriv(double x, double c);

//todo: rename R0 parameter, as its confusing. I think that R0 here is the radius of the tube!
/* Compute initial pressure, density and temperature conditions based on Newton-Raphson's method.  */
ChapmanJougetDetonationSolution SolveChapmanJougetDetonationProblem(const double T0, const double P0, const double ETA, const double S0, const double R, const double GAMMA, const double L_TUBE, const double R0, const double convergenceThreshold = 1.0E-10);



/* End of file */