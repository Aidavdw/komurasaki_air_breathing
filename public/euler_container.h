#pragma once
#include "AuxFunctions.h"

// Describes the properties for a single entity of flow; mostly implemented so that not every flow variable has to be passed through separately in every function.
struct EulerContinuity
{
    EulerContinuity() :
        mass(0), momentumX(0), momentumY(0), energy(0)
    {}

    EulerContinuity(const double mass, const double momentumX, const double momentumY, const double energy) :
        mass(mass),
        momentumX(momentumX),
        momentumY(momentumY),
        energy(energy)
    {}
    
    double mass; // Conservation of mass, first term of the euler continuity.
    double momentumX; // Conservation of momentum in the x-direction, second term of the euler continuity.
    double momentumY; // Conservation of momentum in the y-direction, third term of the euler continuity.
    double energy;  // Conservation of energy, fourth term in the euler continuity.

    inline EulerContinuity operator+ (const EulerContinuity& other) const
    {
        return {mass + other.mass, momentumX + other.momentumX, momentumY + other.momentumY, energy + other.energy};
    }

    inline EulerContinuity operator- (const EulerContinuity& other) const
    {
        return operator+(other*-1);
    }

    inline EulerContinuity operator* (const double scale) const
    {
        return {mass * scale, momentumX * scale, momentumY * scale, energy * scale};
    }
    
    inline EulerContinuity operator/ (const double scale) const
    {
        return operator*(1/scale);
    }

    bool operator== (const EulerContinuity& other) const;

};

inline bool EulerContinuity::operator==(const EulerContinuity& other) const
{
    return (IsCloseToZero(mass - other.mass) && IsCloseToZero(momentumX - other.momentumX) && IsCloseToZero(momentumY - other.momentumY) && IsCloseToZero(energy - other.energy));
}
