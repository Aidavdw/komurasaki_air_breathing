/*
Depending on the case name, initial conditions are specified for all the domains and fields of the total fluid domain. If new cases are defined by the user and new case names are defined, the user should also define an initial condition for the new "reference" case he created. This initial condition should be defined in this function.
*/

const int ARRAY_LEN=14;
const int ARRAY_WIDTH=250; 

void apply_initial_conditions(char *sim_case, int ndom, int *nx, int *ny, double ***x, double ***rho, double ***u, double ***v, double ***p, double ***E, double ***T, double ***H, int nghost, double m_msd, double p1, double u1, double rho1, double m1, double p2, double rho2, double l_exp, double l_tube, double t0, double p0, double m0, double r, double gamma, double *INI_u[ARRAY_WIDTH], double *INI_v[ARRAY_WIDTH], double *INI_rho[ARRAY_WIDTH], double *INI_p[ARRAY_WIDTH])
{
    /* First case defined by user... */
    if (strcmp(sim_case,"det_tube")==0)
    {
        for (int k = 0; k < ndom; ++k)
        {
            double x_center;
            switch(k)
            {
                case 0:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    x_center = 0.5*(x[k][i][nghost]+x[k][i+1][nghost]);
                    if (x_center <= l_exp)
                    {
                        for (int j = nghost; j < ny[k]-nghost; ++j)
                        {
                            p[k][i][j]=p2;
                            rho[k][i][j]=rho2;
                            T[k][i][j]=p2/(rho2*r);
                            u[k][i][j]=0.0;
                            v[k][i][j]=0.0;

                            E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                            H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        }
                    }
                    else if (x_center > l_exp)
                    {
                        for (int j = nghost; j < ny[k]-nghost; ++j)
                        {
                            p[k][i][j]=p2*pow(1.0 - (gamma - 1.0)/(gamma + 1.0)*(1.0-x_center/l_exp),2.0*gamma/(gamma-1.0));
                            rho[k][i][j]=rho2*pow(p[k][i][j]/p2,1.0/gamma);
                            T[k][i][j]=p2/(rho2*r)*pow(rho[k][i][j]/rho2,gamma-1.0);
                            u[k][i][j]=-2.0/(gamma-1.0)*(sqrt(gamma*p2/rho2)-sqrt(gamma*r*T[k][i][j]));
                            v[k][i][j]=0.0;

                            E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                            H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        }
                    }
                }
                break;

                default:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        T[k][i][j]=T0;
                        p[k][i][j]=P0;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                        E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    }
                }
                break;
            }
        }
    }

    /* Detonation tube for the Microwave configuration. This appears to be identical to det_tube. */ 
	else if (strcmp(sim_case,"mw_tube")==0)
	{
		for (int k = 0; k < ndom; ++k)
		{
			double x_center;
			switch(k)
			{
				case 0:
				for (int i = nghost; i < nx[k]-nghost; ++i)
				{
					x_center = 0.5*(x[k][i][nghost]+x[k][i+1][nghost]);
					if (x_center <= l_exp)
					{
						for (int j = nghost; j < ny[k]-nghost; ++j)
						{
	                        p[k][i][j]=p2;
	                        rho[k][i][j]=rho2;
	                        T[k][i][j]=p2/(rho2*r);
	                        u[k][i][j]=0.0;
	                        v[k][i][j]=0.0;

	                        E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
	                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
						}
					}
					else if (x_center > l_exp)
					{
						for (int j = nghost; j < ny[k]-nghost; ++j)
						{
	                        p[k][i][j]=p2*pow(1.0 - (gamma - 1.0)/(gamma + 1.0)*(1.0-x_center/l_exp),2.0*gamma/(gamma-1.0));
	                        rho[k][i][j]=rho2*pow(p[k][i][j]/p2,1.0/gamma);
	                        T[k][i][j]=p2/(rho2*r)*pow(rho[k][i][j]/rho2,gamma-1.0);
	                        u[k][i][j]=-2.0/(gamma-1.0)*(sqrt(gamma*p2/rho2)-sqrt(gamma*r*T[k][i][j]));
	                        v[k][i][j]=0.0;

	                        E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
	                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
						}
					}
				}
				break;

				default:
				for (int i = nghost; i < nx[k]-nghost; ++i)
				{
					for (int j = nghost; j < ny[k]-nghost; ++j)
					{
                        T[k][i][j]=T0;
                        p[k][i][j]=P0;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                        E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
					}
				}
				break;
			}
		}
	}

    /* Microwave Rocket with separated plenum */
	else if (strcmp(sim_case,"plenum_rocket")==0)
	{
        for (int k = 0; k < ndom; ++k)
        {
            double x_center;
            switch(k)
            {
                case 0:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    x_center = 0.5*(x[k][i][nghost]+x[k][i+1][nghost]);
                    if (x_center <= l_exp)
                    {
                        for (int j = nghost; j < ny[k]-nghost; ++j)
                        {
                            p[k][i][j]=p2;
                            rho[k][i][j]=rho2;
                            T[k][i][j]=p2/(rho2*r);
                            u[k][i][j]=0.0;
                            v[k][i][j]=0.0;

                            E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                            H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        }
                    }
                    else if (x_center > l_exp)
                    {
                        for (int j = nghost; j < ny[k]-nghost; ++j)
                        {
                            p[k][i][j]=p2*pow(1.0 - (gamma - 1.0)/(gamma + 1.0)*(1.0-x_center/l_exp),2.0*gamma/(gamma-1.0));
                            rho[k][i][j]=rho2*pow(p[k][i][j]/p2,1.0/gamma);
                            T[k][i][j]=p2/(rho2*r)*pow(rho[k][i][j]/rho2,gamma-1.0);
                            u[k][i][j]=-2.0/(gamma-1.0)*(sqrt(gamma*p2/rho2)-sqrt(gamma*r*T[k][i][j]));
                            v[k][i][j]=0.0;

                            E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                            H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        }
                    }
                }
                break;
                // case 0 ... 2:
                // for (int i = nghost; i < nx[k]-nghost; ++i)
                // {
                //     for (int j = nghost; j < ny[k]-nghost; ++j)
                //     {
                //         T[k][i][j]=T0;
                //         p[k][i][j]=P0*0.4;
                //         u[k][i][j]=0.0;
                //         v[k][i][j]=0.0;

                //         rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                //         E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                //         H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                //     }
                // }
                // break;

                default:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        T[k][i][j]=T0;
                        p[k][i][j]=P0;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                        E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    }
                }
                break;
            }
        }
	}

    /* Modify this condition to obtain convergence of the plenum region */
	else if (strcmp(sim_case,"sup_plenum")==0)
	{
        char *var[] = {"rho","u","v","p","T","E","H"};

		for (int k = 0; k < ndom; ++k)
		{
			switch(k)
			{
				// default:
				// for (int i = nghost; i < nx[k]-nghost; ++i)
				// {
				// 	for (int j = nghost; j < ny[k]-nghost; ++j)
				// 	{
    //                     T[k][i][j]=t0;
    //                     p[k][i][j]=p0;
    //                     u[k][i][j]=m0*sqrt(gamma*r*t0);
    //                     v[k][i][j]=0.0;

    //                     rho[k][i][j]=p[k][i][j]/(T[k][i][j]*r);
    //                     E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
    //                     H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
				// 	}
				// }
				// break;

				// case 1:
    //             printf("Stagnation conditions: P = %f Pa, T = %f K.\n",sup_stag_p(m0,p0,gamma),sup_stag_t(m0,t0,gamma));
				// for (int i = nghost; i < nx[k]-nghost; ++i)
				// {
				// 	for (int j = nghost; j < ny[k]-nghost; ++j)
				// 	{
    //                     T[k][i][j]=sup_stag_t(m0,t0,gamma);
    //                     p[k][i][j]=sup_stag_p(m0,p0,gamma);
    //                     u[k][i][j]=0.0;
    //                     v[k][i][j]=0.0;

    //                     rho[k][i][j]=p[k][i][j]/(T[k][i][j]*r);
    //                     E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
    //                     H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
				// 	}
				// }
				// break;

                default:
                {
                    import_matrix(k,var[0],nx[k],ny[k],nghost,rho[k]);
                    import_matrix(k,var[1],nx[k],ny[k],nghost,u[k]);
                    import_matrix(k,var[2],nx[k],ny[k],nghost,v[k]);
                    import_matrix(k,var[3],nx[k],ny[k],nghost,p[k]);
                    import_matrix(k,var[4],nx[k],ny[k],nghost,T[k]);
                    import_matrix(k,var[5],nx[k],ny[k],nghost,E[k]);
                    import_matrix(k,var[6],nx[k],ny[k],nghost,H[k]);
                    // for (int i = nghost; i < nx[k]-nghost; ++i)
                    // {
                    //     for (int j = nghost; j < ny[k]-nghost; ++j)
                    //     {
                    //         T[k][i][j]=p[k][i][j]/(rho[k][i][j]*r);
                    //         E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                    //         H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    //     }
                    // }
                break;
                }
			}
		}
	}

    /* Microwave Rocket with pre-run plenum region */
    else if (strcmp(sim_case,"sup_plen_rocket")==0)
    {
        for (int k = ndom-1; k > -1; --k)
        {
            double x_center;
            char *var[] = {"rho","u","v","p","E","T","H"};

            switch(k)
            {
                case 0:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    x_center = 0.5*(x[k][i][nghost]+x[k][i+1][nghost]);
                    if (x_center <= l_exp)
                    {
                        for (int j = nghost; j < ny[k]-nghost; ++j)
                        {
                            p[k][i][j]=p2;
                            rho[k][i][j]=rho2;
                            T[k][i][j]=p2/(rho2*r);
                            u[k][i][j]=0.0;
                            v[k][i][j]=0.0;

                            E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                            H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        }
                    }
                    else if (x_center > l_exp)
                    {
                        for (int j = nghost; j < ny[k]-nghost; ++j)
                        {
                            p[k][i][j]=p2*pow(1.0 - (gamma - 1.0)/(gamma + 1.0)*(1.0-x_center/l_exp),2.0*gamma/(gamma-1.0));
                            rho[k][i][j]=rho2*pow(p[k][i][j]/p2,1.0/gamma);
                            T[k][i][j]=p2/(rho2*r)*pow(rho[k][i][j]/rho2,gamma-1.0);
                            u[k][i][j]=-2.0/(gamma-1.0)*(sqrt(gamma*p2/rho2)-sqrt(gamma*r*T[k][i][j]));
                            v[k][i][j]=0.0;

                            E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                            H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        }
                    }
                }
                break;

                case 1 ... 2:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        T[k][i][j]=t0;
                        p[k][i][j]=p0;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*r);
                        E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    }
                }
                break;

                case 3: case 4: case 5: case 6: case 7:
                for (int i = 0; i < 7; ++i)
                {
                    if (i==0)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,rho[k]);
                    }
                    else if (i==1)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,u[k]);
                    }
                    else if (i==2)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,v[k]);
                    }
                    else if (i==3)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,p[k]);
                    }
                    else if (i==4)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,E[k]);
                    }
                    else if (i==5)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,T[k]);
                    }
                    else if (i==6)
                    {
                        import_matrix(k-3,var[i],nx[k],ny[k],nghost,H[k]);
                    }
                }
                // for (int i = nghost; i < nx[k]-nghost; ++i)
                // {
                //     for (int j = nghost; j < ny[k]-nghost; ++j)
                //     {
                //         T[k][i][j]=p[k][i][j]/(rho[k][i][j]*r);
                //         E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                //         H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                //     }
                // }
                break;
            }
        }
    }

    /* Microwave Rocket at ground */
    else if (strcmp(sim_case,"ground_rocket")==0)
    {
        

        for (int k = 0; k < ndom; ++k)
        {
            double x_center;
            switch(k)
            {
                case 0: //the tube case which can be modified with new field obtained from code/measurement
                
                
                for (int i = nghost; i < nx[k]-nghost; ++i)
                { 
                    
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        printf("\nInitial conditions applied...\n");
                        printf("%d", i);
                        printf("%d", j);
                        p[k][i][j]=INI_p[i][j]; 
                        rho[k][i][j]=INI_rho[i][j];
                        T[k][i][j]=INI_p[i][j]/(INI_rho[i][j]*r);
                        u[k][i][j]=INI_u[i][j];
                        v[k][i][j]=INI_v[i][j];

                        E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                        printf("\nInitial conditions applied...\n");
                    }
                    
                  
                }
                break;



                // for (int i = nghost; i < nx[k]-nghost; ++i)
                // {
                //     x_center = 0.5*(x[k][i][nghost]+x[k][i+1][nghost]);
                //     if (x_center <= l_exp)
                //     {
                //         for (int j = nghost; j < ny[k]-nghost; ++j)
                //         {
                //             p[k][i][j]=p2;
                //             rho[k][i][j]=rho2;
                //             T[k][i][j]=p2/(rho2*r);
                //             u[k][i][j]=0.0;
                //             v[k][i][j]=0.0;

                //             E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                //             H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                //         }
                //     }
                //     else if (x_center > l_exp)
                //     {
                //         for (int j = nghost; j < ny[k]-nghost; ++j)
                //         {
                //             p[k][i][j]=p2*pow(1.0 - (gamma - 1.0)/(gamma + 1.0)*(1.0-x_center/l_exp),2.0*gamma/(gamma-1.0));
                //             rho[k][i][j]=rho2*pow(p[k][i][j]/p2,1.0/gamma);
                //             T[k][i][j]=p2/(rho2*r)*pow(rho[k][i][j]/rho2,gamma-1.0);
                //             u[k][i][j]=-2.0/(gamma-1.0)*(sqrt(gamma*p2/rho2)-sqrt(gamma*r*T[k][i][j]));
                //             v[k][i][j]=0.0;

                //             E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                //             H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                //         }
                //     }
                // }
                // break;
                // case 0 :
                // for (int i = nghost; i < nx[k]-nghost; ++i)
                // {
                //     for (int j = nghost; j < ny[k]-nghost; ++j)
                //     {
                //         T[k][i][j]=T0;
                //         p[k][i][j]=P0*0.8;
                //         u[k][i][j]=0.0;
                //         v[k][i][j]=0.0;

                //         rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                //         E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                //         H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                //     }
                // }
                // break;

                default:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        T[k][i][j]=t0;
                        p[k][i][j]=p0;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*r);
                        E[k][i][j]=p[k][i][j]/(gamma-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    }
                }
                break;

                // default:
                // for (int i = nghost; i < nx[k]-nghost; ++i)
                // {
                //     for (int j = nghost; j < ny[k]-nghost; ++j)
                //     {
                //         T[k][i][j]=T0;
                //         p[k][i][j]=P0;
                //         u[k][i][j]=0.0;
                //         v[k][i][j]=0.0;

                //         rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                //         E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                //         H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                //     }
                // }
                // break;
            }
        }
    }
    else if (strcmp(sim_case,"valve_tank")==0)
    {
        for (int k = 0; k < ndom; ++k)
        {
            switch(k)
            {
                case 0:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        T[k][i][j]=T0;
                        p[k][i][j]=P0*0.65;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                        E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    }
                }
                break;

                default:
                for (int i = nghost; i < nx[k]-nghost; ++i)
                {
                    for (int j = nghost; j < ny[k]-nghost; ++j)
                    {
                        T[k][i][j]=T0;
                        p[k][i][j]=P0;
                        u[k][i][j]=0.0;
                        v[k][i][j]=0.0;

                        rho[k][i][j]=p[k][i][j]/(T[k][i][j]*R);
                        E[k][i][j]=p[k][i][j]/(GAMMA-1.0) + 0.5*rho[k][i][j]*(pow(u[k][i][j],2)+pow(v[k][i][j],2));
                        H[k][i][j]=(E[k][i][j]+p[k][i][j])/rho[k][i][j];
                    }
                }
                break;
            }
        }
    }
    else
    {
        fprintf(stderr,"\nERROR: No initial condition matching current case!\n"); exit(0);
    }
}

/* Used for Newton-Raphson method */
double mesh_function(double dx,double alpha,double l,int n)
{
    return dx*(1.0-pow(alpha,n))-l*(1.0-alpha);
}

/* Used for Newton-Raphson method */
double mesh_derivative(double dx,double alpha,double l,int n)
{
    return l-dx*pow(alpha,n-1)*n;
}

/* Calculation of mesh ratios to obtain graded meshes (computationally accurate?) */
double mesh_ratio(int n_cell,double dx0,double l,double grid_ratio)
{
    const double static THRESHOLD = 1.0E-12;

    double alpha_current,alpha_next;
    if (grid_ratio>1.0)
    {
        alpha_next = 1.1;
    }
    else if (grid_ratio<1.0)
    {
        alpha_next = 0.9;
    }
    else if (grid_ratio == 1.0)
    {
        return 1.0;
    }
    int i = 0;

    while(i==0 || fabs(alpha_current-alpha_next) > THRESHOLD)
    {
        alpha_current = alpha_next;
        alpha_next = alpha_current - mesh_function(dx0,alpha_current,l,n_cell)/mesh_derivative(dx0,alpha_current,l,n_cell);
        i++;
    }

    return alpha_next;
}

/* Generate the mesh for the whole domain, more precisely every subdomain of the total fluid domain */
void compute_mesh(int ndom,int *nx, int *ny, double ***x, double ***y, double ***xc, double ***yc, double *xstart, double *ystart, double *xlength, double *ylength, double *xratio, double *yratio, int nghost)
{
	double dx0,dy0,alpha_x,alpha_y;

    for (int k = 0; k < ndom; ++k)
    {
        if (xratio[k]!=1.0 || yratio[k]!=1.0)
        {
            dx0 = xlength[0]/(nx[0]-2*nghost);
            dy0 = ylength[0]/(ny[0]-2*nghost);
        }
        else
        {
            dx0 = xlength[k]/(nx[k]-2*nghost);
            dy0 = ylength[k]/(ny[k]-2*nghost);
        }
        alpha_x = mesh_ratio(nx[k]-2*nghost,dx0,xlength[k],xratio[k]);
        alpha_y = mesh_ratio(ny[k]-2*nghost,dy0,ylength[k],yratio[k]);

        if (xstart[k] >= 0.0 || xratio[k] == 1.0)
        {
            for (int i = 0; i < nx[k]+1; ++i)
            {
                for (int j = 0; j < ny[k]+1; ++j)
                {
                    if (i<=nghost)
                    {
                        x[k][i][j] = xstart[k]+(i-nghost)*dx0;
                    }
                    else if (i > nghost && i < nx[k]-nghost+1)
                    {
                        x[k][i][j] = x[k][i-1][j] + (x[k][i-1][j] - x[k][i-2][j])*alpha_x;
                    }
                    else
                    {
                        x[k][i][j] = 2.0*x[k][i-1][j] - x[k][i-2][j];
                    }
                }
            }
        }
        else if (xstart[k] < 0.0 && xratio[k]!=1.0)
        {
            for (int i = 0; i < nx[k]+1; ++i)
            {
                for (int j = 0; j < ny[k]+1; ++j)
                {
                    int i2 = nx[k]-i;
                    if (i2>=nx[k]-nghost+1)
                    {
                        x[k][i2][j] = xstart[k]+xlength[k]+(i2-nx[k]+nghost)*dx0;
                    }
                    else if (i2 > nghost && i2 < nx[k]-nghost+1)
                    {
                        x[k][i2][j] = x[k][i2+1][j] - (x[k][i2+2][j] - x[k][i2+1][j])*alpha_x;
                    }
                    else
                    {
                        x[k][i2][j] = 2.0*x[k][i2+1][j] - x[k][i2+2][j];
                    }
                }
            }
        }

        for (int i = 0; i < nx[k]+1; ++i)
        {
            for (int j = 0; j < ny[k]+1; ++j)
            {
                if (j<=nghost)
                {
                    y[k][i][j] = ystart[k]+(j-nghost)*dy0;
                }
                else if (j > nghost && j < ny[k]-nghost+1)
                {
                    y[k][i][j] = y[k][i][j-1] + (y[k][i][j-1] - y[k][i][j-2])*alpha_y;
                }
                else
                {
                    y[k][i][j] = 2.0*y[k][i][j-1] - y[k][i][j-2];
                }

                if (i>0 && j>0)
                {
                    // Cell centres
                    xc[k][i-1][j-1] = (x[k][i-1][j-1]+x[k][i][j])/2.0;
                    yc[k][i-1][j-1] = (y[k][i-1][j-1]+y[k][i][j])/2.0;
                }
            }
        }
    }
}

/* END OF FILE */