/*  Flux Vector Splitting schemes used during computation are listed in this file. AUSM-DV formulation is based on the joined use of AUSM-D and AUSM-V, and uses Hanel's dissipative scheme in directions orthogonal to that where a shock is detected (in order to avoid carbuncle phenomenon). An entropy fix is also employed so as to respect the entropy condition and avoid non-physical solutions.

Listed in these files are the following:
    - AUSM: AUSM-DV scheme with entropy fix
    - HANEL: Hanel's dissipative scheme.

Choice between these two schemes is performed directly in the main routine.
*/

/* AUSM-DV scheme at a certain face (R is right, L is left of this face) --> Might be problematic */
void AUSM_DV(double flux[4], const char horOrVer, const double rhoL, const double rhoR, const double uL, const double uR, const double vL, const double vR, const double pL, const double pR, const double HL, const double HR,double r, const double gamma, const double AUSM_K, const double ENTRO_FIX_C)
{
    double s=0;
    double rhou2MV=0,rhou2MD=0,rhouhalf=0,uhalf=0,phalf;
    double uLplus=0,uRminus=0,pLplus=0,pRminus=0;

    double cL=sqrt(gamma*pL/rhoL);
    double cR=sqrt(gamma*pR/rhoR);
    double alphaL=pL/rhoL;
    double alphaR=pR/rhoR;
    double ALPHA_L=2.0*alphaL/(alphaL+alphaR);
    double ALPHA_R=2.0*alphaR/(alphaR+alphaL);
    double cM=fmax(cL,cR);

    // Definition of UL(plus), PL(plus), UR(minus) and PR(minus)
    if (fabs(uL)<=cM)
    {
        uLplus=0.5*(uL + fabs(uL)) + ALPHA_L*(0.25/cM*pow(uL + cM,2) - 0.5*(uL + fabs(uL)));
        pLplus=0.25*pL*pow(uL/cM+1.0,2)*(2.0 - uL/cM);
    }
    else
    {
        uLplus=0.5*(uL+fabs(uL));
        pLplus=0.5*pL*(1.0+fabs(uL)/uL);
    }
    if (fabs(uR)<=cM)
    {
        uRminus=0.5*(uR - fabs(uR)) + ALPHA_R*(-0.25/cM*pow(uR - cM,2) - 0.5*(uR - fabs(uR)));
        pRminus=0.25*pR*pow(uR/cM - 1.0,2)*(2.0+uR/cM);
    }
    else
    {
        uRminus=0.5*(uR - fabs(uR));
        pRminus=0.5*pR*(1.0 - fabs(uR)/uR);
    }

    // Bias function for pressure gradient
    s = 0.5*fmin(1.0,AUSM_K*fabs(pR - pL)/fmin(pR,pL));

    // AUSM-V and AUSM-D momentum flux terms
    uhalf=uLplus+uRminus;
    phalf=pLplus+pRminus;
    rhouhalf=0.5*(uhalf*(rhoR+rhoL) - fabs(uhalf)*(rhoR - rhoL));
    rhou2MD=0.5*(rhouhalf*(uL+uR) - fabs(rhouhalf)*(uR - uL));
    rhou2MV=uLplus*rhoL*uL + uRminus*rhoR*uR;

    // AUSM-DV otherwise
    flux[0] = uLplus*rhoL + uRminus*rhoR;
    flux[1] = (0.5+s)*rhou2MV+(0.5 - s)*rhou2MD+phalf;
    flux[2] = 0.5*(rhouhalf*(vL+vR) - fabs(rhouhalf)*(vR - vL));
    flux[3] = 0.5*(rhouhalf*(HL+HR) - fabs(rhouhalf)*(HR - HL));

    // Entropy fix (numerical dissipation for single expansion waves)
    // Identical for both cases.
    if((uL - cL<0.0 && uR - cR>0.0)&(!(uL+cL<0.0 && uR+cR>0.0)))
    {
        flux[0]=flux[0] - ENTRO_FIX_C*(uR - cR - uL+cL)*(rhoR - rhoL);
        flux[1]=flux[1] - ENTRO_FIX_C*(uR - cR - uL+cL)*(rhoR*uR - rhoL*uL);
        flux[2]=flux[2] - ENTRO_FIX_C*(uR - cR - uL+cL)*(rhoR*vR - rhoL*vL);
        flux[3]=flux[3] - ENTRO_FIX_C*(uR - cR - uL+cL)*(rhoR*HR - rhoL*HL);
    }
    else if((!(uL - cL<0.0 && uR - cR>0.0))&(uL+cL<0.0 && uR+cR>0.0))
    {
        flux[0]=flux[0] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR - rhoL);
        flux[1]=flux[1] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR*uR - rhoL*uL);
        flux[2]=flux[2] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR*vR - rhoL*vL);
        flux[3]=flux[3] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR*HR - rhoL*HL);
    }

    if (horOrVer=='V')
    {
        // Inverse flux components in case of vertical flux
        double tempVar=0;
        tempVar = flux[2];
        flux[2] = flux[1]; // U -> V, V -> -U by 90° rotation
        flux[1] = tempVar;
    }
}

/* Hanel's dissipative scheme. To be used orthogonally to directions where shocks have been detected */
void HANEL(double flux[4], const char horOrVer, const double rhoL, const double rhoR, const double uL, const double uR, const double vL, const double vR, const double pL, const double pR, const double HL, const double HR, const double r, const double gamma, const double ENTRO_FIX_C)
{
    double phalf;
    double uLplus;
    double uRminus;
    double pLplus;
    double pRminus;

    double cL=sqrt(gamma*pL/rhoL); // Speed of sound on the left face.
    double cR=sqrt(gamma*pR/rhoR); // Speed of sound on the right face.
    double alphaL=pL/rhoL;
    double alphaR=pR/rhoR;
    double ALPHA_L=2*alphaL/(alphaL+alphaR);
    double ALPHA_R=2*alphaR/(alphaR+alphaL);
    double cM=fmax(cL,cR);

    // Definition of UL(plus), PL(plus), UR(minus) and PR(minus)
    if (fabs(uL)<=cM)
    {
        uLplus=0.5*(uL+fabs(uL))+ALPHA_L*(0.25/cM*pow(uL+cM,2) - 0.5*(uL+fabs(uL)));
        pLplus=0.25*pL*pow(uL/cM+1,2)*(2 - uL/cM);
    }
    else
    {
        uLplus=0.5*(uL+fabs(uL));
        pLplus=0.5*pL*(1+fabs(uL)/uL);
    }

    if (fabs(uR)<=cM)
    {
        uRminus=0.5*(uR - fabs(uR))+ALPHA_R*( - 0.25/cM*pow(uR - cM,2) - 0.5*(uR - fabs(uR)));
        pRminus=0.25*pR*pow(uR/cM - 1,2)*(2+uR/cM);
    }
    else
    {
        uRminus=0.5*(uR - fabs(uR));
        pRminus=0.5*pR*(1 - fabs(uR)/uR);
    }

    // Hänel scheme for sonic points (shock fix)
    phalf=pLplus+pRminus;
    flux[0] = uLplus*rhoL + uRminus*rhoR;
    flux[1] = uLplus*rhoL*uL + uRminus*rhoR*uR + phalf;
    flux[2] = uLplus*rhoL*vL + uRminus*rhoR*vR;
    flux[3] = uLplus*rhoL*HL + uRminus*rhoR*HR;

    // Entropy fix (numerical dissipation for single expansion waves)
    // Identical for both cases.
    if((uL - cL<0 && uR - cR>0)&(!(uL+cL<0 && uR+cR>0)))
    {
        flux[0]=flux[0] - ENTRO_FIX_C*(uR-cR-uL+cL)*(rhoR - rhoL);
        flux[1]=flux[1] - ENTRO_FIX_C*(uR-cR-uL+cL)*(rhoR*uR - rhoL*uL);
        flux[2]=flux[2] - ENTRO_FIX_C*(uR-cR-uL+cL)*(rhoR*vR - rhoL*vL);
        flux[3]=flux[3] - ENTRO_FIX_C*(uR-cR-uL+cL)*(rhoR*HR - rhoL*HL);
    }
    else if((!(uL-cL<0 && uR-cR>0))&(uL+cL<0 && uR+cR>0))
    {
        flux[0]=flux[0]-ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR - rhoL);
        flux[1]=flux[1] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR*uR - rhoL*uL);
        flux[2]=flux[2] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR*vR - rhoL*vL);
        flux[3]=flux[3] - ENTRO_FIX_C*(uR+cR - uL - cL)*(rhoR*HR - rhoL*HL);
    }

    if (horOrVer=='V')
    {
        // Inverse flux components in case of vertical flux
        double tempVar;
        tempVar = flux[2];
        flux[2] = flux[1]; // U -> V, V -> -U by 90° rotation
        flux[1] = tempVar;
    }
}

/* End of file */