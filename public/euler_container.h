#pragma once
#include "AuxFunctions.h"

// Describes the properties for a single entity of flow; mostly implemented so that not every flow variable has to be passed through separately in every function.
struct EulerContinuity
{
    EulerContinuity() :
        density(0), u(0), v(0), e(0), h(0), t(0)
    {}

    EulerContinuity(const double density, const double xAxisVelocity, const double yAxisVelocity, const double specificInternalEnergy, const double specificEnthalpy, const double staticPressure, const double staticTemperature) :
        density(density),
        u(xAxisVelocity),
        v(yAxisVelocity),
        e(specificInternalEnergy),
        p(staticPressure),
        h(specificEnthalpy),
        t(staticTemperature)
    {}

    double density;
    double u;           // x-axis velocity
    double v;           // y-axis velocity
    double e;           // specific internal energy
    double p;           // static pressure
    double h;           // enthalpy.
    double t;           // static temperature

    inline EulerContinuity operator+ (const EulerContinuity& other) const
    {
        return {density + other.density, u + other.u, v + other.v, e + other.e, h + other.h, p + other.p, t + other.t};
    }

    inline EulerContinuity operator- (const EulerContinuity& other) const
    {
        return operator+(other*-1);
    }

    inline EulerContinuity operator* (const double scale) const
    {
        return {density * scale, u * scale, v * scale, e * scale, h * scale, p*scale, t*scale};
    }
    
    inline EulerContinuity operator/ (const double scale) const
    {
        return operator*(1/scale);
    }

    bool operator== (const EulerContinuity& other) const;

};

inline bool EulerContinuity::operator==(const EulerContinuity& other) const
{
    return (IsCloseToZero(density - other.density) && IsCloseToZero(u - other.u) && IsCloseToZero(v - other.v) && IsCloseToZero(e - other.e) && IsCloseToZero(p - other.p) && IsCloseToZero(h - other.h) && IsCloseToZero(t - other.t));
}
