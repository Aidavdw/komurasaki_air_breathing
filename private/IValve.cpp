#include "IValve.h"

#include <cassert>
#include <stdexcept>

#include "domain.h"

IValve::~IValve()
{
}

IValve::IValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, double positionAlongBoundary) :
    boundary_(boundary),
    positionAlongBoundary_(positionAlongBoundary),
    outOfDomain_(outOfDomain),
    intoDomain_(intoDomain)
{
    pos = intoDomain->PositionAlongBoundaryToCoordinate(boundary, positionAlongBoundary, 0);
}

void IValve::CalculateFlow()
{
#ifdef _DEBUG
    assert(sourceCellsIndices_.size() == sourceTermBuffer_.size());
    assert(intoDomain_ != nullptr);
    if (sourceCellsIndices_.empty())
        throw std::logic_error("No source cells are set for this valve.");

    for (const auto& sourceTerms : sourceTermBuffer_)
    {
        if (!(sourceTerms == EulerContinuity()))
            throw std::logic_error("Source term is non-empty."); 
    }
#endif
    
    FillSourceTermBuffer();
}

void IValve::OnRegister()
{
    throw std::logic_error("OnRegister() is not overridden for this type of valve!");
}


void IValve::UpdateValveState()
{
    throw std::logic_error("Update() is not overridden for this type of valve!");
}

void IValve::SetInitialConditions()
{
    throw std::logic_error("SetInitialConditions() is not overridden for this type of valve!");
}

void IValve::FillSourceTermBuffer()
{
    throw std::logic_error("PopulateValveDeltaBuffer() is not overridden for this type of valve!");
}

void IValve::AddCachedTermsToDomainConservationEquations() const
{
#ifdef _DEBUG
    assert(sourceCellsIndices_.size() == sourceTermBuffer_.size());
    assert(sinkCellsIndices_.size() == sinkTermBuffer_.size());
    assert(sinkCellsIndices_.size() == sourceCellsIndices_.size());
    assert(intoDomain_ != nullptr);
    if (sourceCellsIndices_.empty())
        throw std::logic_error("No source cells are set for this valve.");

    // The total sum of all fluxes must be zero
    EulerContinuity totalSource;
    EulerContinuity totalSink;
    for (auto& sourceTermInCell : sourceTermBuffer_)
        totalSource = totalSource + sourceTermInCell;
    for (auto& sinkTermInCell : sinkTermBuffer_)
        totalSink = totalSink + sinkTermInCell;
    if (!((totalSource - totalSink) == EulerContinuity()))
        throw std::logic_error("conservation laws are being broken, as the total source term is not equal to the total sink term for this valve!");

#endif
    
    for (size_t i = 0; i < sourceCellsIndices_.size(); i++)
    {
        // 1. The source cells in the in-domain
        const CellIndex& cixSource = sourceCellsIndices_[i];
        const EulerContinuity sourceTerms = sourceTermBuffer_[i];
        
#ifdef _DEBUG
        assert(intoDomain_->ValidateCellIndex(sourceCellsIndices_.at(i), false));
        if (sourceTerms == EulerContinuity())
            throw std::logic_error("Source term is not initialised.");
#endif


        // todo: confirm that the source term is already divided by the amount of cells. if not, add division by source cells.
        intoDomain_->eulerConservationTerms[0](cixSource) += sourceTerms.mass;
        intoDomain_->eulerConservationTerms[1](cixSource) += sourceTerms.momentumX;
        intoDomain_->eulerConservationTerms[2](cixSource) += sourceTerms.momentumY;
        intoDomain_->eulerConservationTerms[3](cixSource) += sourceTerms.energy;

        // 2. The sink cells in the outOfDomain_
        const CellIndex& cixSink = sinkCellsIndices_[i];
        const EulerContinuity sinkTerms = sinkTermBuffer_[i];
        
#ifdef _DEBUG
        assert(outOfDomain_->ValidateCellIndex(sinkCellsIndices_.at(i), false));
        if (sinkTerms == EulerContinuity())
            throw std::logic_error("sink term is not initialised.");
#endif

        outOfDomain_->eulerConservationTerms[0](cixSource) -= sinkTerms.mass;
        outOfDomain_->eulerConservationTerms[1](cixSource) -= sinkTerms.momentumX;
        outOfDomain_->eulerConservationTerms[2](cixSource) -= sinkTerms.momentumY;
        outOfDomain_->eulerConservationTerms[3](cixSource) -= sinkTerms.energy;
    }
    
}

void IValve::EmptyBuffer()
{
#ifdef _DEBUG
    assert(sourceCellsIndices_.size() == sourceTermBuffer_.size());
    assert(sinkCellsIndices_.size() == sinkTermBuffer_.size());
    assert(sinkCellsIndices_.size() == sourceCellsIndices_.size());
#endif
    
    for (size_t i = 0; i < sourceCellsIndices_.size() ; i++)
    {
        sourceTermBuffer_[i] = EulerContinuity();
        sinkTermBuffer_[i] = EulerContinuity();
    }
}

std::string IValve::ToString() const
{
    return "Valve in '" + intoDomain_->name + "' at " + FaceToString(boundary_) + " dist="+ std::to_string(positionAlongBoundary_);
}

