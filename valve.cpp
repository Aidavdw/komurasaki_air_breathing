#include "valve.h"
#include "domain.h"

Valve::Valve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary) :
    intoDomain(intoDomain),
    boundary(boundary),
    positionAlongBoundary(positionAlongBoundary)
{
    
}

void Valve::OnRegister()
{
    // Accessing the source cells:
    //intoDomain->T.main.FlattenIndexOnBoundary(boundary,sourceStartIndexOnBoundary + idxBelowEndIndex)
}




int Valve::GetSourceTermCount() const
{
    sourceEndIndexOnBoundary - sourceStartIndexOnBoundary;
}
