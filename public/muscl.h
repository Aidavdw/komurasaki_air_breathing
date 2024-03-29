﻿#pragma once
#include "domain_enums.h"

enum class EFluxLimiterType
{
    NONE,
    MIN_MOD,
    SUPER_BEE,
    VAN_ALBADA_ONE,
    VAN_ALBADA_TWO
};

// Represents the values of a variable along its edges, as MUSCL interpolated with its neighbours. [0] is interpolated to the inside 凹, [1] is to the outside 凸.
struct MUSCLBuffer
{
    double left[2] = {0,0};
    double right[2] = {0,0};
    double top[2] = {0,0};
    double bottom[2] = {0,0};

    double GetAt(const EFace face, EAxisDirection axisDirection) const;
    double& At(EFace face, EAxisDirection axisDirection);
};

// Returns the MUSCL interpolated value of centre, given the values of its left neighbour m1, and two right neighbours, p1 and p2.
double MUSCLInterpolate(const double m1, const double centre, const double p1, const double p2, const EAxisDirection sideToInterpolateTo, const double bias, const EFluxLimiterType fluxLimiterType);

double ApplyFluxLimiter(const double r, const EFluxLimiterType fluxLimiterType);

// Typical values of flux limiters
// +----------+---------+----------+-------------+--------------+
// | input r  | minmod  | superbee |   albada1   |   albada2    |
// +----------+---------+----------+-------------+--------------+
// |    -2000 |       0 |        0 | 0.99949975  | -0.001       |
// |       -5 |       0 |        0 | 0.769230769 | -0.384615385 |
// |  -0.0002 |       0 |        0 | -0.00019996 | -0.0004      |
// |        0 |       0 |        0 | 0           | 0            |
// |  0.00002 | 0.00002 |  0.00004 | 2.00004E-05 | 4E-05        |
// |        5 |       1 |        2 | 1.153846154 | 0.384615385  |
// | 99999999 |       1 |        2 | 1.00000001  | 2E-08        |
// +----------+---------+----------+-------------+--------------+


/* Minmod limiter function. r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)) */
double MinMod(const double r);

/* Super-Bee limiter function. r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)) */
double SuperBee(const double r);

/* Van Albada limiter function (version 1: TVD). r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i) */
double VanAlbadaOne(const double r);

// Van Albada limiter function (version 2: not TVD). r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)). */
double VanAlbadaTwo(const double r);