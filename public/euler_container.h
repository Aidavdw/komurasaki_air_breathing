#pragma once

struct EulerContinuity
{
    EulerContinuity() :
        density(0), u(0), v(0), e(0), h(0)
    {}

    EulerContinuity(const double density, const double u, const double v, const double e, const double enthalpy, const double p) :
        density(density),
        u(u),
        v(v),
        e(e),
        p(p),
        h(enthalpy)
    {}

    double density;
    double u;           // x-axis velocity
    double v;           // y-axis velocity
    double e;           // specific internal energy
    double p;           // pressure
    double h;           // enthalpy.

    inline EulerContinuity operator+ (const EulerContinuity& other) const
    {
        return {density + other.density, u + other.u, v + other.v, e + other.e, h + other.h, p + other.p};
    }

    inline EulerContinuity operator- (const EulerContinuity& other) const
    {
        return operator+(other*-1);
    }

    inline EulerContinuity operator* (const double scale) const
    {
        return {density * scale, u * scale, v * scale, e * scale, h * scale, p*scale};
    }
    
    inline EulerContinuity operator/ (const double scale) const
    {
        return operator*(1/scale);
    }
};
