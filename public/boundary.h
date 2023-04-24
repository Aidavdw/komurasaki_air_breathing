#pragma once
#include "domain_enums.h"

struct Domain;

// Describes a boundary that a domain has four of (considering it is 2d, and all domains are rectangular
struct Boundary
{
    Boundary() :
        boundaryType(EBoundaryType::NOT_SET),
        connectedBoundary(nullptr)
    {}

    explicit Boundary(const EBoundaryType boundaryType) :
        boundaryType(boundaryType),
        connectedBoundary(nullptr)
    {}

    // This boundary itself:
    //Domain* domain; // The domain that this boundary is a part of
    EBoundaryType boundaryType;

    // If it's connected to another boundary too:
    Boundary* connectedBoundary;
};
