#pragma once

struct EulerContinuity
{
    EulerContinuity() :
        density(0), u(0), v(0), e(0)
    {}

    EulerContinuity(const double density, const double u, const double v, const double e) :
        density(density),
        u(u),
        v(v),
        e(e)
    {}

    double density;
    double u;           // x-axis velocity
    double v;           // y-axis velocity
    double e;           // specific internal energy

    inline EulerContinuity operator+ (const EulerContinuity& other) const
    {
        return {density + other.density, u + other.u, v + other.v, e + other.e};
    }

    inline EulerContinuity operator- (const EulerContinuity& other) const
    {
        return operator+(other*-1);
    }

    inline EulerContinuity operator* (const double scale) const
    {
        return {density * scale, u * scale, v * scale, e * scale};
    }
    
    inline EulerContinuity operator/ (const double scale) const
    {
        return operator*(1/scale);
    }
};
