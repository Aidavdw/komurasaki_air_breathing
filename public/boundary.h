#pragma once
#include "domain_enums.h"

struct Domain;

// Describes a boundary that a domain has four of (considering it is 2d, and all domains are rectangular
struct Boundary
{
    Boundary() :
        boundaryType(EBoundaryCondition::NOT_SET),
        connectedBoundary(nullptr),
        domain(nullptr)
    {}

    explicit Boundary(const EBoundaryCondition boundaryType) :
        boundaryType(boundaryType),
        connectedBoundary(nullptr),
        domain(nullptr)
    {}

    // This boundary itself:
    //Domain* domain; // The domain that this boundary is a part of
    EBoundaryCondition boundaryType;

    // If it's connected to another boundary too:
    Boundary* connectedBoundary;
    Domain* domain;
};
