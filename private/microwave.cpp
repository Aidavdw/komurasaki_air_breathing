#include "microwave.h"

#include <cassert>

#include "domain.h"
#include <cmath>
#include <stdexcept>
#include "AuxFunctions.h"


double eval_msd_function(double x, double c)
{
    return std::sqrt(c / x + 1) + std::sqrt(c / x) - x;
}

double eval_msd_deriv(double x, double c)
{
    return -c / std::pow(x, 2) * (0.5 * std::sqrt(c / x + 1) + 0.5 * std::sqrt(c / x)) - 1;
}

ChapmanJougetDetonationSolution SolveChapmanJougetDetonationProblem(const double temperatureAmbient, const double pressureAmbient, const double eta, const double s0, const double idealGasConstant, const double gamma, const double lengthOfCombustionTube, const double radiusOfCombustionTube, const double convergenceThreshold)
{

    const static int MAX_ITERS = 200;

    double rho0 = pressureAmbient / (idealGasConstant * temperatureAmbient); // Ambient density
    double a0 = std::sqrt(gamma * idealGasConstant * temperatureAmbient); // Ambient speed of sound
    double powerDensity = eta * s0 / (M_PI * radiusOfCombustionTube * radiusOfCombustionTube); // Mean power density
    
    const double cFactor = (std::pow(gamma, 2) - 1) / (2 * std::pow(a0, 3) * rho0) * powerDensity; // Constant factor for Newton-Raphson method

    double Mcur = 1, Mnext = 0.0;
    Mnext = eval_msd_function(1.0,cFactor);
    int iters_performed = 0;

    while (std::abs(Mcur - Mnext) > convergenceThreshold)
    {
        if (iters_performed > MAX_ITERS)
            throw std::runtime_error("Chapman-Jouget detonation took more than the maximum amount of iterations to solve. Are you sure you're doing this right?");
        Mcur = Mnext;
        Mnext = Mcur - eval_msd_function(Mcur, cFactor) / eval_msd_deriv(Mcur, cFactor);

#ifdef _DEBUG
        assert(!std::isnan(Mnext));
#endif
        iters_performed++;
    }

    // Do the same calculation, but now with ETA factor for comparison
    // A: This appears to not have been implemented fully. Leaving here for the time being
    //TODO: either remove this eta case, or implement it.
    /*
    const double C_eta = (std::pow(GAMMA, 2) - 1) / (2 * std::pow(a0, 3) * rho0) * 0.49 * Sd;
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
    ChapmanJougetDetonationSolution sol;
    sol.m1 = (std::pow(Mnext, 2) - 1.0) / (1.0 + gamma * std::pow(Mnext, 2));

    // Post-detonation conditions, immediately behind the detonation front (subscript 1). These are required for calculating the plateau conditions, and afterwards discarded. The calculation of the values behind the shockwave is x-location dependent, and is defined in ChapmanJougetDetonationSolution::FieldPropertiesAtPosition()
    CellValues postDetonation;           // The values of the flow after the expansion
    postDetonation.p = (1.0 + gamma * std::pow(Mnext, 2)) / (gamma + 1) * pressureAmbient;
    postDetonation.u = a0 * Mnext * (std::pow(Mnext, 2.0) - 1.0) / std::pow(Mnext, 2.0) / (gamma + 1.0);
    postDetonation.density = (1.0 + gamma) * std::pow(Mnext, 2) / (1.0 + gamma * std::pow(Mnext, 2)) * rho0;

    // Post-expansion wave conditions subscript 2
    CellValues postExpansion;           // The values of the flow after the expansion
    postExpansion.p = std::pow(1.0 - 0.5 * (gamma - 1) * (sol.m1), 2.0 * gamma / (gamma - 1)) * (postDetonation.p);
    postExpansion.density = std::pow((1.0 - 0.5 * (gamma - 1) * (sol.m1)), 2.0 / (gamma - 1)) * (postDetonation.density);
    double a2 = (1.0 - 0.5 * (gamma - 1.0) * (sol.m1)) / (gamma + 1.0) * (std::pow(Mnext, 2) * gamma + 1.0) / std::pow(Mnext, 2) * Mnext * a0;
    postExpansion.t = postExpansion.p / (postExpansion.density * idealGasConstant);
    postExpansion.u = 0;
    postExpansion.v = 0;
    postExpansion.e = postExpansion.p / (gamma - 1.0) + 0.5 * postExpansion.density * (std::pow(postExpansion.u, 2) + std::pow(postExpansion.v, 2)); // since u & v are 0, these terms drop out entirely.
    postExpansion.h = (postExpansion.e + postExpansion.p) / postExpansion.density;
    sol.m_msd = Mnext;
    // Computing position of expansion wave front and rear
    sol.l_exp = (lengthOfCombustionTube / Mnext / a0) * a2;
    sol.detonation_velocity = Mnext * a0;

    sol.postDetonation = postDetonation;
    sol.postExpansion = postExpansion;
    return sol;
}

void InitialiseDomainFromChapmanJougetDetonationSolution(Domain* domain, const ChapmanJougetDetonationSolution& sol, const double gamma)
{
    //Todo: Take this parameter out, so that it can also be made standing up.
    const double tubeRadius = domain->size[1]; // For now, assume a laying down tube. 
    // As the solution is 1D, the solution is the same across all the y - coordinates, and is only different for x coordinates. Hence, iterate over the x- coordinates and set it immediately for all y coordinates.
    for (int xIndex = 0; xIndex < domain->amountOfCells[0]; xIndex++)
    {
        double xPos = domain->localCellCenterPositions[0].at(xIndex);
        CellValues cellValues = sol.FieldPropertiesAtPosition(xPos, gamma, tubeRadius);
#ifdef _DEBUG
        // Make sure we're not setting fields to 0.
        assert(!IsCloseToZero(cellValues.density));
        //assert(!IsCloseToZero(cellValues.v)); this will probably be zero lol
        //assert(!IsCloseToZero(cellValues.u)); this will probably be zero lol
        assert(!IsCloseToZero(cellValues.p));
        assert(!IsCloseToZero(cellValues.e));
        assert(!IsCloseToZero(cellValues.h));
            
#endif
        
        for (int yIndex = 0; yIndex < domain->amountOfCells[1]; yIndex++)
        {
            domain->rho.currentTimeStep(xIndex, yIndex) = cellValues.density;
            domain->u.currentTimeStep(xIndex, yIndex) = (cellValues.u);
            domain->v.currentTimeStep(xIndex, yIndex) = (cellValues.v);
            domain->p.currentTimeStep(xIndex, yIndex) = (cellValues.p);
            domain->E.currentTimeStep(xIndex, yIndex) = (cellValues.e);
            domain->T.currentTimeStep(xIndex, yIndex) = (cellValues.t);
            domain->H.currentTimeStep(xIndex, yIndex) = (cellValues.h);
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
        c.p = postExpansion.p* std::pow(1.0 - (gamma - 1.0) / (gamma + 1.0) * (1.0 - xPosition / l_exp), 2.0 * gamma / (gamma - 1.0));
        c.density = postExpansion.density * std::pow(postExpansion.p / c.p, 1.0 / gamma);
        c.t = postExpansion.p / (postExpansion.density * tubeRadius) * std::pow(c.density / postExpansion.density, gamma - 1.0);
        c.u = -2.0 / (gamma - 1.0) * (sqrt(gamma * postExpansion.p / postExpansion.density) - std::sqrt(gamma * tubeRadius * c.t));
        c.v = 0.0;

        c.e = c.p / (gamma - 1.0) + 0.5 * c.density * (std::pow(c.u, 2) + std::pow(c.v, 2));
        c.h= (c.e + c.p) / c.density;

        return c;
    }
}
