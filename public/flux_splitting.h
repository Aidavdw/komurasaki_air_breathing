#pragma once
/* This file contains the two functions for flux splitting used in solving the euler equations.
 */
#include "euler_container.h"

// forward declarations
class Domain;

// Flux splitting scheme used if the flow is (super)sonic.
EulerContinuity HanelFluxSplitting(const EulerContinuity& l, const EulerContinuity& r, const double gamma, const double entropyFix);  

// Flux splitting scheme used if the flow is sub-sonic.
EulerContinuity AUSMDVFluxSplitting(const EulerContinuity& leftContinuity, const EulerContinuity& rightContinuity, const double gamma, const double AUSMkFactor, const double entropyFix);