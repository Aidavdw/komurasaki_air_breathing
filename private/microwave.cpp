#include "microwave.h"
#include "domain.h"
#include <math.h>
#include <stdexcept>


double eval_msd_function(double x, double c)
{
    return sqrt(c / x + 1) + sqrt(c / x) - x;
}

double eval_msd_deriv(double x, double c)
{
    return -c / pow(x, 2) * (0.5 * sqrt(c / x + 1) + 0.5 * sqrt(c / x)) - 1;
}

ChapmanJougetDetonationSolution SolveChapmanJougetDetonationProblem(const double T0, const double P0, const double ETA, const double S0, const double R_ideal, const double GAMMA, const double L_TUBE, const double RadiusOfCombustionTube, const double convergenceThreshold)
{
    CellValues postExpansion;           // The values of the flow after the expansion

    // todo: This value was defined in main.c, not sure where it came from. Move to different scope?
    const static int M_PI = 4.0 * atan(1.0);
    const static int MAX_ITERS = 200;

    double rho0 = P0 / (R_ideal * T0); // Ambient density
    double a0 = sqrt(GAMMA * R_ideal * T0); // Ambient speed of sound
    double Sd = ETA * S0 / (M_PI * RadiusOfCombustionTube * RadiusOfCombustionTube); // Mean power density

    //printf("Microwave surface power density is: %f MW/cm^2.\n",Sd/1.0E6/100.0/100.0);
    // printf("Ionization velocity is: %f m/s.\n",4190.0*Sd/1.0E6/100.0/100.0-14.9);
    // Constant factor for Newton-Raphson method
    const double C = (pow(GAMMA, 2) - 1) / (2 * pow(a0, 3) * rho0) * Sd;

    double Mcur = 1, Mnext = 0.0;
    //Mnext = eval_msd_function(1.0,C);
    int iters_performed = 0;

    while (fabs(Mcur - Mnext) > convergenceThreshold)
    {
        if (iters_performed > MAX_ITERS)
            throw std::runtime_error("Chapman-Jouget detonation took more than the maximum amount of iterations to solve. Are you sure you're doing this right?");
        Mcur = Mnext;
        Mnext = Mcur - eval_msd_function(Mcur, C) / eval_msd_deriv(Mcur, C);
        iters_performed++;
    }



    // Do the same calculation, but now with ETA factor for comparison
    // A: This appears to not have been implemented fully. Leaving here for the time being
    //TODO: either remove this eta case, or implement it.
    /*
    const double C_eta = (pow(GAMMA, 2) - 1) / (2 * pow(a0, 3) * rho0) * 0.49 * Sd;
    Mcur = 0.0;
    double Meta = 0.0;
    Meta = eval_msd_function(1.0, C_eta);
    while (fabs(Mcur - Meta) > convergenceThreshold)
    {
        Mcur = Meta;
        Meta = Mcur - eval_msd_function(Mcur, C_eta) / eval_msd_deriv(Mcur, C_eta);
    }
    //printf("Detonation velocity with ETA=0.49 is: %f m/s.\n",Meta*a0);
    */

    // This is equation 3.11 from Florian (2017)
    double m1 = (pow(Mnext, 2) - 1.0) / (1.0 + GAMMA * pow(Mnext, 2));

    // Post-detonation conditions, immediately behind the detonation front (subscript 1). These are required for calculating the plateau conditions, and afterwards discarded. The calculation of the values behind the shockwave is x-location dependent, and is defined in ChapmanJougetDetonationSolution::FieldPropertiesAtPosition()
    double p1 = (1.0 + GAMMA * pow(Mnext, 2)) / (GAMMA + 1) * P0;
    double u1 = a0 * Mnext * (pow(Mnext, 2.0) - 1.0) / pow(Mnext, 2.0) / (GAMMA + 1.0);
    double rho1 = (1.0 + GAMMA) * pow(Mnext, 2) / (1.0 + GAMMA * pow(Mnext, 2)) * rho0;


    // Post-expansion wave conditions subscript 2
    postExpansion.p = pow(1.0 - 0.5 * (GAMMA - 1) * (m1), 2.0 * GAMMA / (GAMMA - 1)) * (p1);
    postExpansion.rho = pow((1.0 - 0.5 * (GAMMA - 1) * (m1)), 2.0 / (GAMMA - 1)) * (rho1);
    double a2 = (1.0 - 0.5 * (GAMMA - 1.0) * (m1)) / (GAMMA + 1.0) * (pow(Mnext, 2) * GAMMA + 1.0) / pow(Mnext, 2) * Mnext * a0;
    postExpansion.T = postExpansion.p / (postExpansion.rho * R_ideal);
    postExpansion.u = 0;
    postExpansion.v = 0;
    postExpansion.E = postExpansion.p / (GAMMA - 1.0) + 0.5 * postExpansion.rho * (pow(postExpansion.u, 2) + pow(postExpansion.v, 2)); // since u & v are 0, these terms drop out entirely.
    postExpansion.H = (postExpansion.E + postExpansion.p) / postExpansion.rho;



    ChapmanJougetDetonationSolution solution = ChapmanJougetDetonationSolution(postExpansion);
    solution.m_msd = Mnext;
    solution.m1 = m1;
    // Computing position of expansion wave front and rear
    solution.l_exp = (L_TUBE / Mnext / a0) * a2;
    solution.detonation_velocity = Mnext * a0;

    // printf("\nMean Microwave Power is: %f W/m^2.\n",Sd);
    //printf("Detonation velocity is: %f m/s (M_MSD = %f).\n", Mnext * a0, Mnext);
    //printf("Ambient conditions are: \nP0 = %f Pa, \nRHO0 = %f, \na0 = %f.", P0, rho0, a0);
    //printf("Post-MSD conditions are: \nM1 = %f, \nP1 = %f Pa, \nU1 = %f m/s, \nRHO1 = %f.\n", *m1, *p1, *u1, *rho1);
    //printf("Post-expansion conditions are: \nP2 = %f Pa, \nRHO2 = %f, \na2 = %f.\n", *p2, *rho2, a2);
    //printf("Tube length is %f m. The tail of the expansion region is at %f m.\n", L_TUBE, *l_exp);

    return solution;
}

void InitialiseDomainFromChapmanJougetDetonationSolution(Domain* domain, const ChapmanJougetDetonationSolution& sol, const double gamma)
{
    //Todo: Take this parameter out, so that it can also be made standing up.
    const double tubeRadius = domain->size[1]; // For now, assume a laying down tube. 
    // As the solution is 1D, the solution is the same across all the y - coordinates, and is only different for x coordinates. Hence, iterate over the x- coordinates and set it immediately for all y coordinates.
    for (int xIndex = 0; xIndex < domain->size[0]; xIndex++)
    {
        double xPos = domain->localCellCenterPositions[0].GetAt(xIndex,0);
        CellValues cellValues = sol.FieldPropertiesAtPosition(xPos, gamma, tubeRadius);
        
        for (int yIndex = 0; yIndex < domain->size[1]; yIndex++)
        {
            domain->rho.SetAllToValue(cellValues.rho);
            domain->u.SetAllToValue(cellValues.u);
            domain->v.SetAllToValue(cellValues.v);
            domain->p.SetAllToValue(cellValues.p);
            domain->E.SetAllToValue(cellValues.E);
            domain->T.SetAllToValue(cellValues.T);
            domain->E.SetAllToValue(cellValues.E);
        }
    }
}

CellValues ChapmanJougetDetonationSolution::FieldPropertiesAtPosition(const double xPosition, const double gamma, const double tubeRadius) const
{
    // See if this is the plateau region or if this is the region where stuff gradually increases.
    /*     |           /
    *  P,T |          /
    *  U   |         /
    *  rho | -------+
    *      |
    *      +------------
    *        post   ^  inside expansion
    *             l_exp
    */
    
    if (xPosition <= l_exp)
    {
        return postExpansion;
    }
    else // It's in the expansion region- Note that the solution is only valid up until the actual x position at the end of the tube.
    {
        CellValues c;
        c.p = postExpansion.p* pow(1.0 - (gamma - 1.0) / (gamma + 1.0) * (1.0 - xPosition / l_exp), 2.0 * gamma / (gamma - 1.0));
        c.rho = postExpansion.rho * pow(postExpansion.p / c.p, 1.0 / gamma);
        c.T = postExpansion.p / (postExpansion.rho * tubeRadius) * pow(c.rho / postExpansion.rho, gamma - 1.0);
        c.u = -2.0 / (gamma - 1.0) * (sqrt(gamma * postExpansion.p / postExpansion.rho) - sqrt(gamma * tubeRadius * c.T));
        c.v = 0.0;

        c.E = c.p / (gamma - 1.0) + 0.5 * c.rho * (pow(c.u, 2) + pow(c.v, 2));
        c.H= (c.E + c.p) / c.rho;

        return c;
    }
}
