#include "IValve.h"

#include <cassert>
#include <stdexcept>

#include "domain.h"

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
    
    FillBuffer();
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

void IValve::FillBuffer()
{
    throw std::logic_error("PopulateValveDeltaBuffer() is not overridden for this type of valve!");
}

void IValve::AddBufferTermsToSourceCells(const EFieldQuantityBuffer buffer)
{
#ifdef _DEBUG
    assert(sourceCellsIndices_.size() == sourceTermBuffer_.size());
    assert(intoDomain_ != nullptr);
    if (sourceCellsIndices_.empty())
        throw std::logic_error("No source cells are set for this valve.");
#endif

    TwoDimensionalArray& rhom = intoDomain_->rho.bufferMap.at(buffer);
    TwoDimensionalArray& um = intoDomain_->u.bufferMap.at(buffer);
    TwoDimensionalArray& vm = intoDomain_->v.bufferMap.at(buffer);
    TwoDimensionalArray& pm = intoDomain_->p.bufferMap.at(buffer);
    TwoDimensionalArray& em = intoDomain_->E.bufferMap.at(buffer);
    TwoDimensionalArray& tm = intoDomain_->T.bufferMap.at(buffer);
    TwoDimensionalArray& hm = intoDomain_->H.bufferMap.at(buffer);
    
    for (size_t i = 0; i < sourceCellsIndices_.size(); i++)
    {
        const CellIndex& cix = sourceCellsIndices_[i];
        const EulerContinuity sourceTerms = sourceTermBuffer_[i];
        
#ifdef _DEBUG
        assert(intoDomain_->ValidateCellIndex(sourceCellsIndices_.at(i), false));
        if (sourceTerms == EulerContinuity())
            throw std::logic_error("Source term is not initialised.");
#endif

        rhom(cix) = sourceTerms.density;
        um(cix) = sourceTerms.u;
        vm(cix) = sourceTerms.v;
        pm(cix) = sourceTerms.p;
        em(cix) = sourceTerms.e;
        hm(cix) = sourceTerms.h;

        //todo: set temperature based on the other values of state
        tm(cix) = 0;
    }
    
}

void IValve::EmptyBuffer()
{
#ifdef _DEBUG
    assert(sourceCellsIndices_.size() == sourceTermBuffer_.size());
#endif
    
    for (size_t i = 0; i < sourceCellsIndices_.size() ; i++)
        sourceTermBuffer_[i] = EulerContinuity();
}

