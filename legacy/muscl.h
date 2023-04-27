
/*  Here are functions that are used in the framework of MUSCL interpolation. Piecewise parabolic MUSCL method is used for theorical 3rd order space accuracy in continuous regions.
*/

/* Minmod limiter function. r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)) */
double minmod(double r)
{
    return fmax(0.0,fmin(1.0,r));
}

/* Super-Bee limiter function. r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)) */
double superbee(double r)
{
    return fmax(fmax(0.0,fmin(2*r,1.0)),fmax(0.0,fmin(r,2.0)));
}

/* Van Albada limiter function (version 1: TVD). r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i) */
double van_albada1(double r)
{
    return (pow(r,2)+r)/(pow(r,2)+1);
}

// Van Albada limiter function (version 2: not TVD). r(i) is defined as (u(i)-u(i-1))/(u(i+1)-u(i)). */
double van_albada2(double r)
{
    return 2.0*r/(pow(r,2)+1);
}

/* Calculate the limited flux using specified limiter function */
double fluxLimiter(double r, int limiterName)
{
    switch(limiterName)
    {
        case -1: return 1.0;
        case 0 : return minmod(r);
        case 1 : return superbee(r);
        case 2 : return van_albada1(r);
        case 3 : return van_albada2(r);
        default : printf("\nERROR: Invalid flux-limiter function: '%d'!",limiterName); exit(0);break;
    }
}

/* Calculate MUSCL interpolation based on points (-1,0,+1,+2) */
double MUSCL(const double m1, const double centre, const double p1, const double p2, const char RorL, const double bias, const int limiterName)
{
    double var=0;
    double deltaplus,deltamin;

    switch(RorL)
    {
        case 'R' :  deltaplus = p2 - p1;
                    deltamin = p1 - centre;
                    break;
        case 'L' :  deltaplus = p1 - centre;
                    deltamin = centre - m1;
                    break;
        default :   deltaplus = 0.0;
                    deltamin = 0.0;
                    printf("\nERROR: Specified MUSCL keyword is not 'left' or 'right': %c",RorL); break; exit(0);
    }

    // Check for possible zero division in smooth fields
    double ratio;
    double inverse;
    double eps = 1.0E-7;
    ratio = (deltamin*deltaplus+eps)/(deltaplus*deltaplus+eps);
    inverse = (deltamin*deltaplus+eps)/(deltamin*deltamin+eps);


    switch(RorL)
    {
        case 'R' :  var = p1 -0.25*((1-bias)*deltaplus*fluxLimiter(ratio,limiterName)+(1+bias)*deltamin*fluxLimiter(inverse,limiterName));
                    break;
        case 'L' :  var = centre + 0.25*((1-bias)*deltamin*fluxLimiter(inverse,limiterName)+(1+bias)*deltaplus*fluxLimiter(ratio,limiterName));
                    break;
        default :   var=0.0;
                    printf("\nERROR: Specified MUSCL keyword is not 'left' or 'right': %c",RorL); break; exit(0);
    }

    return var;
}

/* End of file */