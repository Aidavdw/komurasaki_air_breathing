#include "parameters.h"

/*  This file contains functions that are related with Microwave Beaming and MSD theory:
    - eval_msd_function: Evaluates the value of function to solve to derive MSD Mach number.
    - eval_msd_deriv: Evaluates the value of the derivative of the former function.
    - solve_msd: Estimation of MSD Mach number based on Newton Raphson method. Returns the number of iterations required for sufficient convergence.

*/

struct ChapmanJougetDetonationSolution
{
    double m_msd = 0;
    double p1 = 0;
    double u1 = 0;
    double rho = 0;
    double m1 = 0;
    double p2 = 0;
    double rho2 = 0;
    double l_exp = 0;
};

/* Function to solve to find MSD Mach number */
double eval_msd_function(double x, double c);

/* Derivative of function to solve to find MSD Mach number */
double eval_msd_deriv(double x, double c);


// T0,P0,ETA,S0,R,GAMMA,L_TUBE,R0
/* Compute initial pressure, density and temperature conditions based on Newton-Raphson's method */
int solve_MSD(double* m_msd, double* p1, double* u1, double* rho1, double* m1, double* p2, double* rho2, double* l_exp);



/* End of file */