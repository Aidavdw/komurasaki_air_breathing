#include "valve.h"
#include "domain.h"
#include <stdexcept>

Valve::Valve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary) :
    intoDomain(intoDomain),
    boundary(boundary),
    positionAlongBoundary(positionAlongBoundary)
{
    
}

void Valve::OnRegister()
{
    throw std::logic_error("OnRegister() is not overridden for this type of valve!");
    // Accessing the source cells:
    //intoDomain->T.main.FlattenIndexOnBoundary(boundary,sourceStartIndexOnBoundary + idxBelowEndIndex)
}

void Valve::SetSourceTerms()
{
    throw std::logic_error("SetSourceTerms() is not overridden for this type of valve!");
}
void Valve::GetAveragePressure() const
{
    throw std::logic_error("GetAveragePressure() is not overridden for this type of valve!");
}

