#include "flux_splitting.h"

#include <algorithm>
#include <cmath>

#include "AuxFunctions.h"
#include "cell_values_container.h"

EulerContinuity HanelFluxSplitting(CellValues l, CellValues r, const bool bIsVertical, const double gamma, const double entropyFix)
{
#ifdef _DEBUG
    l.Validate();
    r.Validate();
#endif

    // Just trust me on this one. This is how Florian did it apparently.
    if (bIsVertical)
    {
        SwapValues(l.u, l.v);
        SwapValues(r.u, r.v);
    }
    
    const double alphaLPartial=l.p/l.density;
    const double alphaRPartial=r.p/r.density;
    const double alphaL=2*alphaLPartial/(alphaLPartial+alphaRPartial);
    const double alphaR=2*alphaRPartial/(alphaRPartial+alphaLPartial);
    
    const double cLeft = std::sqrt(gamma * l.p / l.density);
    const double cRight = std::sqrt(gamma * r.p / r.density);
    const double maxSpeedOfSound =  std::max(cLeft, cRight);

    double uLplus;
    double pLplus;
    if (l.u <= maxSpeedOfSound)
    {
        uLplus=0.5*(l.u+abs(l.u))+alphaL*(0.25/maxSpeedOfSound*std::pow(l.u+maxSpeedOfSound,2)-0.5*(l.u+abs(l.u)));
        pLplus=0.25*l.p*std::pow(l.u/maxSpeedOfSound+1,2)*(2-l.u/maxSpeedOfSound);
    }
    else
    {
        uLplus=0.5*(l.u+abs(l.u));
        pLplus=0.5*l.p*(1+abs(l.u)/l.u);
    }

    double uRminus;
    double pRminus;
    if (abs(r.u)<=maxSpeedOfSound)
    {
        uRminus=0.5*(r.u-abs(r.u))+alphaR*(-0.25/maxSpeedOfSound*std::pow(r.u-maxSpeedOfSound,2)-0.5*(r.u-abs(r.u)));
        pRminus=0.25*r.p*std::pow(r.u/maxSpeedOfSound-1,2)*(2+r.u/maxSpeedOfSound);
    }
    else
    {
        uRminus=0.5*(r.u-abs(r.u));
        pRminus=0.5*r.p*(1-abs(r.u)/r.u);
    }

    // Hänel scheme for sonic points (shock fix)
    const double pHalf=pLplus+pRminus;
    EulerContinuity flux;
    flux.mass = uLplus*l.density+uRminus*r.density;
    flux.momentumX = uLplus*l.density*l.u+uRminus*r.density*r.u+pHalf;
    flux.momentumY = uLplus*l.density*l.v+uRminus*r.density*r.v;
    flux.energy = uLplus*l.density*l.h+uRminus*r.density*r.h;

    // Entropy fix (numerical dissipation for single expansion waves)
    if ( (l.u-cLeft<0 && r.u-cRight>0) && !(l.u+cLeft<0 && r.u+cRight>0))
    {
        flux.mass = flux.mass - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density - l.density);
        flux.momentumX = flux.momentumX - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density * r.u - l.density * l.u);
        flux.momentumY = flux.momentumY - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density * r.v - l.density * l.v);
        flux.energy = flux.energy - entropyFix * (r.u - cRight - l.u + cLeft) * (r.density * r.h - l.density * l.h);
    }
    else if( !(l.u-cLeft<0 && r.u-cRight>0) && (l.u+cLeft<0 && r.u+cRight>0) )
    {
        flux.mass = flux.mass - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density - l.density);
        flux.momentumX = flux.momentumX - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density * r.u - l.density * l.u);
        flux.momentumY = flux.momentumY - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density * r.v - l.density * l.v);
        flux.energy = flux.energy - entropyFix * (r.u + cRight - l.u - cLeft) * (r.density * r.h - l.density * l.h);
    }
    
    // If it's a vertical flux, then the u and v are in the local reference frame, and hence they must be inverted to get it in the r-y reference frame.
    if (bIsVertical)
        SwapValues(flux.momentumX, flux.momentumY);
    return flux;
}

EulerContinuity AUSMDVFluxSplitting(CellValues l, CellValues r, const bool bIsVertical, const double gamma, const double AUSMkFactor, const double entropyFix)
{
#ifdef _DEBUG
    l.Validate();
    r.Validate();
#endif

    // Just trust me on this one. This is how Florian did it apparently.
    if (bIsVertical)
    {
        SwapValues(l.u, l.v);
        SwapValues(r.u, r.v);
    }
    
    const double alphaLPartial=l.p/l.density;
    const double alphaRPartial=r.p/r.density;
    const double alphaL=2*alphaLPartial/(alphaLPartial+alphaRPartial);
    const double alphaR=2*alphaRPartial/(alphaRPartial+alphaLPartial);
    
    const double cL = std::sqrt(gamma * l.p / l.density);
    const double cR = std::sqrt(gamma * r.p / r.density);
    const double maxSpeedOfSound =  std::max(cL, cR);

    double uLplus,pLplus;
    // Definition of UL(plus), PL(plus), UR(minus) and PR(minus)
    if (abs(l.u)<=maxSpeedOfSound)
    {
        uLplus=0.5*(l.u+abs(l.u))+alphaL*(0.25/maxSpeedOfSound*std::pow(l.u+maxSpeedOfSound,2)-0.5*(l.u+abs(l.u)));
        pLplus=0.25*l.p*std::pow(l.u/maxSpeedOfSound+1.0,2)*(2.0-l.u/maxSpeedOfSound);
    }
    else
    {
        uLplus=0.5*(l.u+abs(l.u));
        pLplus=0.5*l.p*(1.0+abs(l.u)/l.u);
    }

    double uRminus,pRminus;
    if (abs(r.u)<=maxSpeedOfSound)
    {
        uRminus=0.5*(r.u-abs(r.u))+alphaR*(-0.25/maxSpeedOfSound*std::pow(r.u-maxSpeedOfSound,2)-0.5*(r.u-abs(r.u)));
        pRminus=0.25*r.p*std::pow(r.u/maxSpeedOfSound-1.0,2)*(2.0+r.u/maxSpeedOfSound);
    }
    else
    {
        uRminus=0.5*(r.u-abs(r.u));
        pRminus=0.5*r.p*(1.0-abs(r.u)/r.u);
    }

    // Bias function for pressure gradient
    double s = 0.5 * std::min(1.0, AUSMkFactor * abs(r.p - l.p) / std::min(r.p, l.p));

    // AUSM-V and AUSM-D momentum flux terms
    double uhalf = uLplus + uRminus;
    double phalf = pLplus + pRminus;
    const double rhoUHalf = 0.5 * (uhalf * (r.density + l.density) - abs(uhalf) * (r.density - l.density));
    const double rhoU2Md = 0.5 * (rhoUHalf * (l.u + r.u) - abs(rhoUHalf) * (r.u - l.u));
    const double rhou2Mv = uLplus * l.density * l.u + uRminus * r.density * r.u;

    EulerContinuity flux;
    // AUSM-DV otherwise
    flux.mass = uLplus*l.density+uRminus*r.density;
    flux.momentumX = (0.5+s)*rhou2Mv+(0.5-s)*rhoU2Md+phalf;
    flux.momentumY = 0.5*(rhoUHalf*(l.v+r.v)-abs(rhoUHalf)*(r.v-l.v));
    flux.energy = 0.5*(rhoUHalf*(l.h+r.h)-abs(rhoUHalf)*(r.h-l.h));

    // Entropy fix (numerical dissipation for single expansion waves)
    if((l.u-cL<0.0 && r.u-cR>0.0)&(!(l.u+cL<0.0 && r.u+cR>0.0)))
    {
        flux.mass = flux.mass - entropyFix * (r.u - cR - l.u + cL) * (r.density - l.density);
        flux.momentumX = flux.momentumX - entropyFix * (r.u - cR - l.u + cL) * (r.density * r.u - l.density * l.u);
        flux.momentumY = flux.momentumY - entropyFix * (r.u - cR - l.u + cL) * (r.density * r.v - l.density * l.v);
        flux.energy = flux.energy - entropyFix * (r.u - cR - l.u + cL) * (r.density * r.h - l.density * l.h);
    }
    else if((!(l.u-cL<0.0 && r.u-cR>0.0))&(l.u+cL<0.0 && r.u+cR>0.0))
    {
        flux.mass = flux.mass - entropyFix * (r.u + cR - l.u - cL) * (r.density - l.density);
        flux.momentumX = flux.momentumX - entropyFix * (r.u + cR - l.u - cL) * (r.density * r.u - l.density * l.u);
        flux.momentumY = flux.momentumY - entropyFix * (r.u + cR - l.u - cL) * (r.density * r.v - l.density * l.v);
        flux.energy = flux.energy - entropyFix * (r.u + cR - l.u - cL) * (r.density * r.h - l.density * l.h);
    }

    // If it's a vertical flux, then the u and v are in the local reference frame, and hence they must be inverted to get it in the r-y reference frame.
    if (bIsVertical)
        SwapValues(flux.momentumX, flux.momentumY);

    return flux;
}
