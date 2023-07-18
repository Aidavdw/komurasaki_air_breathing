#include "muscl.h"

#include <stdexcept>
#include <cmath>

#include "AuxFunctions.h"

#define MAX_FACTOR_DEVIATION_BEFORE_THROWING_ERROR 1.1

double MUSCLBuffer::GetAt(const EFace face, EAxisDirection axisDirection) const
{
    const bool i = (axisDirection == EAxisDirection::POSITIVE);
    switch (face)
    {
    case LEFT: return left[i];
    case RIGHT: return right[i];
    case TOP: return top[i];
    case BOTTOM: return bottom[i];
    default:
        throw std::logic_error("MUSCL buffer getting not implement for this face direction.");
    }
}

double& MUSCLBuffer::At(EFace face, EAxisDirection axisDirection)
{
    const bool i = (axisDirection == EAxisDirection::POSITIVE);
    switch (face)
    {
        case LEFT: return left[i];
        case RIGHT: return right[i];
        case TOP: return top[i];
        case BOTTOM: return bottom[i];
        default:
            throw std::logic_error("MUSCL buffer getting not implement for this face direction.");
    }
}

double MUSCLInterpolate(const double m1, const double centre, const double p1, const double p2, const EAxisDirection sideToInterpolateTo, const double bias, const EFluxLimiterType fluxLimiterType)
{
    //todo: investigate how this works if the cells are not equally sized; how does accumulation work in the ghost cells?

    /*
     *  Interpolates the values at center using the values that are more around it.
     *  Depending on whether left or right is chosen, the entire calculation is the same, but shifted by an index. The final answer is also negated.
     *  Arguably the calculation could be simplified, but that would be at the cost of legibility.
     *
     *  |                           X        
     *  |   X               X                 
     *  |       L   X   R                      
     *  |                                     
     *  |                                    
     *  +---+-------+-------+-------+-------  
     *     m1    center     p1      p2          
     *  
     */

    
    double solution=0;
    double deltaPlus;   // The difference of the two given future values the positive side, aka the side to interpolate to.
    double deltaMin;    // The difference in the negative side, aka relative to the cell behind it.

    switch(sideToInterpolateTo)
    {
    case EAxisDirection::POSITIVE :
        deltaPlus = p2 - p1;
        deltaMin = p1 - centre;
        break;
    case EAxisDirection::NEGATIVE :
        deltaPlus = p1 - centre;
        deltaMin = centre - m1;
        break;
    default :
        throw std::logic_error("MUSCL interpolation is not implemented for this type of side");
    }

    // In a smooth field, it is possible that either deltaPlus or DeltaMin = 0. Dividing by 0 is impossible, so just add a tiny term to it so that the result is finite. This will always be picked up by the flux limiter, and it will be maxed.
    double machineEpsilon = 1.0E-7;
    const double ratioPre = (deltaMin * deltaPlus + machineEpsilon) / (deltaPlus * deltaPlus + machineEpsilon);
    const double inversePre = (deltaMin * deltaPlus + machineEpsilon) / (deltaMin * deltaMin + machineEpsilon);

    switch(sideToInterpolateTo)
    {
    case EAxisDirection::POSITIVE:
        solution = p1 -0.25*((1-bias)*deltaPlus*ApplyFluxLimiter(ratioPre,fluxLimiterType)+(1+bias)*deltaMin*ApplyFluxLimiter(inversePre,fluxLimiterType));
        break;
    case EAxisDirection::NEGATIVE:
        solution = centre + 0.25*((1-bias)*deltaMin*ApplyFluxLimiter(inversePre,fluxLimiterType)+(1+bias)*deltaPlus*ApplyFluxLimiter(ratioPre,fluxLimiterType));
        break;
    default :
        throw std::logic_error("MUSCL interpolation is not implemented for this type of side");
    }

    
#ifdef _DEBUG
    // It's interpolation, so it must be somewhere in between all the values.
    const double max = std::max({m1, centre, p1, p2}) + machineEpsilon;
    double maxLimit = (max > 0) ? max * MAX_FACTOR_DEVIATION_BEFORE_THROWING_ERROR : max / MAX_FACTOR_DEVIATION_BEFORE_THROWING_ERROR;
    if (solution > maxLimit)
        throw std::logic_error("MUSCL interpolation should yield a value between all the values, not above.");
    const double min = std::min({m1, centre, p1, p2}) - machineEpsilon;
    double minLimit = (min > 0) ? min / MAX_FACTOR_DEVIATION_BEFORE_THROWING_ERROR : min * MAX_FACTOR_DEVIATION_BEFORE_THROWING_ERROR;
    if (solution < minLimit)
        throw std::logic_error("MUSCL interpolation should yield a value between all the values, not above.");
#endif
    

    return solution;
    
}

double ApplyFluxLimiter(const double r, const EFluxLimiterType fluxLimiterType)
{
    switch (fluxLimiterType)
    {
    case EFluxLimiterType::NONE: return r;
    case EFluxLimiterType::MIN_MOD: return MinMod(r);
    case EFluxLimiterType::SUPER_BEE: return SuperBee(r);
    case EFluxLimiterType::VAN_ALBADA_ONE: return VanAlbadaOne(r);
    case EFluxLimiterType::VAN_ALBADA_TWO: return VanAlbadaTwo(r);
    default:
        throw std::logic_error("Flux limiter typeis not implemented!");
    }
}

/* Minmod limiter function. r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)) */
double MinMod(const double r)
{
    return fmax(0.0,fmin(1.0,r));
}

/* Super-Bee limiter function. r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)) */
double SuperBee(const double r)
{
    return fmax(fmax(0.0,fmin(2*r,1.0)),fmax(0.0,fmin(r,2.0)));
}

/* Van Albada limiter function (version 1: TVD). r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i) */
double VanAlbadaOne(const double r)
{
    return (std::pow(r,2)+r)/(std::pow(r,2)+1);
}

// Van Albada limiter function (version 2: not TVD). r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)). */
double VanAlbadaTwo(const double r)
{
    return 2.0*r/(std::pow(r,2)+1);
}