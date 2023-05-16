#include "IValve.h"
#include <stdexcept>

IValve::IValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, double positionAlongBoundary) :
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

double IValve::GetMassFlowRate() const
{
    throw std::logic_error("GetMassFlowRate() is not overridden for this type of valve!");
}

void IValve::SetInitialConditions()
{
    throw std::logic_error("SetInitialConditions() is not overridden for this type of valve!");
}

void IValve::PopulateValveDeltaBuffer()
{
    throw std::logic_error("ApplySourceToDomain() is not overridden for this type of valve!");
}

