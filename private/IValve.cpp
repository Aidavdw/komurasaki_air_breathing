#include "IValve.h"
#include <stdexcept>

IValve::IValve(Domain* intoDomain, Domain* outOfDomain, const EBoundaryLocation boundary, double positionAlongBoundary) :
    boundary_(boundary),
    positionAlongBoundary_(positionAlongBoundary),
    outOfDomain_(outOfDomain),
    intoDomain_(intoDomain)
{
    
}

void IValve::OnRegister()
{
    throw std::logic_error("OnRegister() is not overridden for this type of valve!");
}


void IValve::Update()
{
    throw std::logic_error("Update() is not overridden for this type of valve!");
}

void IValve::GetMassFlowRate() const
{
    throw std::logic_error("GetMassFlowRate() is not overridden for this type of valve!");
}

void IValve::SetInitialConditions()
{
    throw std::logic_error("SetInitialConditions() is not overridden for this type of valve!");
}

