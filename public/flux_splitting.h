#pragma once
#include "euler_container.h"

class Domain;

EulerContinuity HanelFluxSplitting(const EulerContinuity& l, const EulerContinuity& r, const double gamma, const double AUSMkFactor, const double entropyFix);

EulerContinuity AUSMDVFluxSplitting(const EulerContinuity& leftContinuity, const EulerContinuity& rightContinuity, const double gasConstant, const double gamma, const double AUSMkFactor, const double entropyFix);