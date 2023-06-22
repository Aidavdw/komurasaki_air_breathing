#pragma once
#include "beam_section.h"

// Describes the physical geometry of a reed valve. Made into a little struct to make it easier to pass through.
struct ReedValveGeometry
{
    EBeamProfile beamProfile = EBeamProfile::STRAIGHT_DOUBLE_TAPERED;
    double freeLength = 0;
    double rootThickness = 0;
    double tipThickness = 0;
    double rootWidth = 0;
    double tipWidth = 0;
};

struct ReedValveEmpiricalParameters
{
    double naturalFrequency = 0;
    double rayleighDampingAlpha = 0;    // Alpha coef. for Rayleigh damping (alpha*M + beta*K)
    double rayleighDampingBeta = 0;     // Beta coef. for Rayleigh damping (alpha*M + beta*K)

    // todo: include flow damping factors c1,c2,c3. They're buried somewhere in the code...
    double dampingC1 = 0;
    double dampingC2 = 0;
    double dampingC3 = 0;

    double holeFactor = 0.9;
};

struct MaterialProperties
{
    double youngsModulus = 0;
    double density = 0;
};