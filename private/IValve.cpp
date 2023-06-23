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

void IValve::AddCachedTermsToSourceCells(const EFieldQuantityBuffer bufferToWriteInto)
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

    // References to the domains it will be added into
    TwoDimensionalArray& rhom = intoDomain_->rho.Buffer(bufferToWriteInto);
    TwoDimensionalArray& um = intoDomain_->u.Buffer(bufferToWriteInto);
    TwoDimensionalArray& vm = intoDomain_->v.Buffer(bufferToWriteInto);
    TwoDimensionalArray& pm = intoDomain_->p.Buffer(bufferToWriteInto);
    TwoDimensionalArray& em = intoDomain_->E.Buffer(bufferToWriteInto);
    TwoDimensionalArray& tm = intoDomain_->T.Buffer(bufferToWriteInto);
    TwoDimensionalArray& hm = intoDomain_->H.Buffer(bufferToWriteInto);

    // References to the domains it will be subtracted from
    TwoDimensionalArray& rhos = outOfDomain_->rho.Buffer(bufferToWriteInto);
    TwoDimensionalArray& us = outOfDomain_->u.Buffer(bufferToWriteInto);
    TwoDimensionalArray& vs = outOfDomain_->v.Buffer(bufferToWriteInto);
    TwoDimensionalArray& ps = outOfDomain_->p.Buffer(bufferToWriteInto);
    TwoDimensionalArray& es = outOfDomain_->E.Buffer(bufferToWriteInto);
    TwoDimensionalArray& ts = outOfDomain_->T.Buffer(bufferToWriteInto);
    TwoDimensionalArray& hs = outOfDomain_->H.Buffer(bufferToWriteInto);
    
    for (size_t i = 0; i < sourceCellsIndices_.size(); i++)
    {
        // The source cells in the in-domain
        const CellIndex& cixSource = sourceCellsIndices_[i];
        const EulerContinuity sourceTerms = sourceTermBuffer_[i];
        
#ifdef _DEBUG
        assert(intoDomain_->ValidateCellIndex(sourceCellsIndices_.at(i), false));
        if (sourceTerms == EulerContinuity())
            throw std::logic_error("Source term is not initialised.");
#endif

        rhom(cixSource) = sourceTerms.density;
        um(cixSource) = sourceTerms.u;
        vm(cixSource) = sourceTerms.v;
        pm(cixSource) = sourceTerms.p;
        em(cixSource) = sourceTerms.e;
        hm(cixSource) = sourceTerms.h;

        //todo: set temperature based on the other values of state
        tm(cixSource) = 0;

        // The sink cells in the outOfDomain_
        const CellIndex& cixSink = sinkCellsIndices_[i];
        const EulerContinuity sinkTerms = sinkTermBuffer_[i];
        
#ifdef _DEBUG
        assert(outOfDomain_->ValidateCellIndex(sinkCellsIndices_.at(i), false));
        if (sinkTerms == EulerContinuity())
            throw std::logic_error("sink term is not initialised.");
#endif

        rhos(cixSink) = sinkTerms.density;
        us(cixSink) = sinkTerms.u;
        vs(cixSink) = sinkTerms.v;
        ps(cixSink) = sinkTerms.p;
        es(cixSink) = sinkTerms.e;
        hs(cixSink) = sinkTerms.h;

        //todo: set temperature based on the other values of state
        ts(cixSink) = 0;
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

