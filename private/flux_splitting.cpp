#include "flux_splitting.h"

#include <algorithm>
#include <cmath>

EulerContinuity HanelFluxSplitting(const EulerContinuity& l, const EulerContinuity& r,
                                   const double gamma, const double AUSMkFactor, const double entropyFix)
{
    const double alphaLPartial=l.p/l.density;
    const double alphaRPartial=r.p/r.density;
    const double alphaL=2*alphaLPartial/(alphaLPartial+alphaRPartial);
    const double alphaR=2*alphaRPartial/(alphaRPartial+alphaLPartial);
    
    const double cLeft = sqrt(gamma * l.p / l.density);
    const double cRight = sqrt(gamma * r.p / r.density);
    const double maxSpeedOfSound =  std::max(cLeft, cRight);

    double uLplus;
    double pLplus;
    if (l.u <= maxSpeedOfSound)
    {
        uLplus=0.5*(l.u+abs(l.u))+alphaL*(0.25/maxSpeedOfSound*pow(l.u+maxSpeedOfSound,2)-0.5*(l.u+abs(l.u)));
        pLplus=0.25*l.p*pow(l.u/maxSpeedOfSound+1,2)*(2-l.u/maxSpeedOfSound);
    }
    else
    {
        uLplus=0.5*(l.u+fabs(l.u));
        pLplus=0.5*l.p*(1+fabs(l.u)/l.u);
    }

    double uRminus;
    double pRminus;
    if (fabs(r.u)<=maxSpeedOfSound)
    {
        uRminus=0.5*(r.u-fabs(r.u))+alphaR*(-0.25/maxSpeedOfSound*pow(r.u-maxSpeedOfSound,2)-0.5*(r.u-fabs(r.u)));
        pRminus=0.25*r.p*pow(r.u/maxSpeedOfSound-1,2)*(2+r.u/maxSpeedOfSound);
    }
    else
    {
        uRminus=0.5*(r.u-fabs(r.u));
        pRminus=0.5*r.p*(1-fabs(r.u)/r.u);
    }

    // Hänel scheme for sonic points (shock fix)
    const double pHalf=pLplus+pRminus;
    EulerContinuity flux;
    flux.density = uLplus*l.density+uRminus*r.density;
    flux.u = uLplus*l.density*l.u+uRminus*r.density*r.u+pHalf;
    flux.v = uLplus*l.density*l.v+uRminus*r.density*r.v;
    flux.e = uLplus*l.density*l.h+uRminus*r.density*r.h;

    // Entropy fix (numerical dissipation for single expansion waves)
    if ( (l.u-cLeft<0 && r.u-cRight>0) && !(l.u+cLeft<0 && r.u+cRight>0))
    {
        flux.density = flux.density - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density - l.density);
        flux.u = flux.u - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density * r.u - l.density * l.u);
        flux.v = flux.v - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density * r.v - l.density * l.v);
        flux.e = flux.e - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density * r.h - l.density * l.h);
    }
    else if( !(l.u-cLeft<0 && r.u-cRight>0) && (l.u+cLeft<0 && r.u+cRight>0) )
    {
        flux.density = flux.density - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density - l.density);
        flux.u = flux.u - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density * r.u - l.density * l.u);
        flux.v = flux.v - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density * r.v - l.density * l.v);
        flux.e = flux.e - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density * r.h - l.density * l.h);
    }

    return flux;
}
