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
    double rayleighDampingAlpha = 0;
    double rayleighDampingBeta = 0;
};