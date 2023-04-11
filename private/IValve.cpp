#include "IValve.h"
#include "domain.h"
#include <stdexcept>

IValve::IValve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary) :
    intoDomain(intoDomain),
    boundary(boundary),
    positionAlongBoundary(positionAlongBoundary)
{
    
}

void IValve::OnRegister()
{
    throw std::logic_error("OnRegister() is not overridden for this type of valve!");
}

void IValve::SetSourceTerms()
{
    throw std::logic_error("SetSourceTerms() is not overridden for this type of valve!");
}
void IValve::GetAveragePressure() const
{
    throw std::logic_error("GetAveragePressure() is not overridden for this type of valve!");
}
void IValve::Update()
{
    throw std::logic_error("Update() is not overridden for this type of valve!");
}

