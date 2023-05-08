#pragma once

struct EulerContainer
{
    EulerContainer(const double density, const double u, const double v, const double e) :
        density(density),
        u(u),
        v(v),
        e(e)
    {}

    double density;
    double u;           // x-axis velocity
    double v;           // y-axis velocity
    double e;           // specific internal energy
};
