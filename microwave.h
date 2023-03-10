#include "parameters.h"

/*  This file contains functions that are related with Microwave Beaming and MSD theory:
    - eval_msd_function: Evaluates the value of function to solve to derive MSD Mach number.
    - eval_msd_deriv: Evaluates the value of the derivative of the former function.
    - solve_msd: Estimation of MSD Mach number based on Newton Raphson method. Returns the number of iterations required for sufficient convergence.

*/

/* Function to solve to find MSD Mach number */
double eval_msd_function(double x, double c)
{
    return sqrt(c/x+1)+sqrt(c/x)-x;
}


/* Derivative of function to solve to find MSD Mach number */
double eval_msd_deriv(double x, double c)
{
    return -c/pow(x,2)*(0.5*sqrt(c/x+1)+0.5*sqrt(c/x))-1;
}

// T0,P0,ETA,S0,R,GAMMA,L_TUBE,R0
/* Compute initial pressure, density and temperature conditions based on Newton-Raphson's method */
int solve_MSD(double *m_msd, double *p1, double *u1, double *rho1, double *m1, double *p2, double *rho2, double *l_exp)
{
    const double static THRESHOLD = 1.0E-10;

    double rho0 = P0/(R*T0); // Ambient density
    double a0 = sqrt(GAMMA*R*T0); // Ambient speed of sound
    double Sd = ETA*S0/(M_PI*R0*R0); // Mean power density

    //printf("Microwave surface power density is: %f MW/cm^2.\n",Sd/1.0E6/100.0/100.0);
    // printf("Ionization velocity is: %f m/s.\n",4190.0*Sd/1.0E6/100.0/100.0-14.9);
    // Constant factor for Newton-Raphson method
    const double C = (pow(GAMMA,2)-1)/(2*pow(a0,3)*rho0)*Sd;

    double Mcur=1, Mnext=0.0;
    //Mnext = eval_msd_function(1.0,C);
    int iters_performed = 0;

    while(fabs(Mcur-Mnext) > THRESHOLD)
    {
        Mcur = Mnext;
        Mnext = Mcur - eval_msd_function(Mcur,C)/eval_msd_deriv(Mcur,C);
        iters_performed++;
    }

    // Return detonation Mach number by reference
    *m_msd = Mnext;

    // Do the same calculation, but now with ETA factor for comparison
    const double C_eta = (pow(GAMMA,2)-1)/(2*pow(a0,3)*rho0)*0.49*Sd;
    Mcur=0.0;
    double Meta = 0.0;
    Meta = eval_msd_function(1.0,C_eta);
    while(fabs(Mcur-Meta) > THRESHOLD)
    {
        Mcur = Meta;
        Meta = Mcur - eval_msd_function(Mcur, C_eta)/eval_msd_deriv(Mcur, C_eta);
    }
    //printf("Detonation velocity with ETA=0.49 is: %f m/s.\n",Meta*a0);

    // Post-detonation conditions
    *m1 = (pow(Mnext,2)-1.0)/(1.0+GAMMA*pow(Mnext,2));
    *p1 = (1.0+GAMMA*pow(Mnext,2))/(GAMMA+1)*P0;
    *u1 = a0*Mnext*(pow(Mnext,2.0)-1.0)/pow(Mnext,2.0)/(GAMMA+1.0);
    *rho1 = (1.0+GAMMA)*pow(Mnext,2)/(1.0+GAMMA*pow(Mnext,2))*rho0;

    // Post-expansion wave conditions
    *p2 = pow(1.0-0.5*(GAMMA-1)*(*m1),2.0*GAMMA/(GAMMA-1))*(*p1);
    *rho2 = pow((1.0-0.5*(GAMMA-1)*(*m1)),2.0/(GAMMA-1))*(*rho1);
    double a2 = (1.0-0.5*(GAMMA-1.0)*(*m1))/(GAMMA+1.0)*(pow(Mnext,2)*GAMMA+1.0)/pow(Mnext,2)*Mnext*a0;

    // Computing position of expansion wave front and rear
    *l_exp = (L_TUBE/Mnext/a0)*a2;


    // printf("\nMean Microwave Power is: %f W/m^2.\n",Sd);
    printf("Detonation velocity is: %f m/s (M_MSD = %f).\n",Mnext*a0,Mnext);
    printf("Ambient conditions are: \nP0 = %f Pa, \nRHO0 = %f, \na0 = %f.",P0,rho0,a0);
    printf("Post-MSD conditions are: \nM1 = %f, \nP1 = %f Pa, \nU1 = %f m/s, \nRHO1 = %f.\n",*m1,*p1,*u1,*rho1);
    printf("Post-expansion conditions are: \nP2 = %f Pa, \nRHO2 = %f, \na2 = %f.\n",*p2,*rho2,a2);
    printf("Tube length is %f m. The tail of the expansion region is at %f m.\n",L_TUBE,*l_exp);
    return iters_performed;
}



/* End of file */