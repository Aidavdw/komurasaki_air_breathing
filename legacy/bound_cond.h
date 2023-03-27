 /* In this file are listed all the boundary conditions used during computation and all related functions. The user can specify new boundary conditions and add them to the existing set (modifications also required in the main file).
*/


/* Connected boundary condition between two fluid domains. This condition ensures that even if the whole domain is defined by blocks, all the blocks act as one consistent domain. */
void connected_condition(int nx1,int nx2,int ny1,int ny2,double **rho1,double **rho2,
    double **p1,double **p2,double **u1,double **u2,double **v1,double **v2,double **H1,
    double **H2,double **E1, double **E2,double **T1,double **T2,char location,int nghost)
{
    switch(location)
    {
        case 'L':
            // Domain 2 is at the left of Domain 1
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < ny1-nghost; ++i)
            {
                rho1[0][i]=rho2[nx2-nghost-2][i];
                u1[0][i]=u2[nx2-nghost-2][i];
                v1[0][i]=v2[nx2-nghost-2][i];
                p1[0][i]=p2[nx2-nghost-2][i];
                H1[0][i]=H2[nx2-nghost-2][i];
                E1[0][i]=E2[nx2-nghost-2][i];
                T1[0][i]=T2[nx2-nghost-2][i];

                rho1[1][i]=rho2[nx2-nghost-1][i];
                u1[1][i]=u2[nx2-nghost-1][i];
                v1[1][i]=v2[nx2-nghost-1][i];
                p1[1][i]=p2[nx2-nghost-1][i];
                H1[1][i]=H2[nx2-nghost-1][i];
                E1[1][i]=E2[nx2-nghost-1][i];
                T1[1][i]=T2[nx2-nghost-1][i];
            } break;

        case 'R':         
            // Domain 2 is at the right of Domain 1
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < ny1-nghost; ++i)
            {
                rho1[nx1-nghost][i]=rho2[nghost][i];
                u1[nx1-nghost][i]=u2[nghost][i];
                v1[nx1-nghost][i]=v2[nghost][i];
                p1[nx1-nghost][i]=p2[nghost][i];
                H1[nx1-nghost][i]=H2[nghost][i];
                E1[nx1-nghost][i]=E2[nghost][i];
                T1[nx1-nghost][i]=T2[nghost][i];

                rho1[nx1-nghost+1][i]=rho2[nghost+1][i];
                u1[nx1-nghost+1][i]=u2[nghost+1][i];
                v1[nx1-nghost+1][i]=v2[nghost+1][i];
                p1[nx1-nghost+1][i]=p2[nghost+1][i];
                H1[nx1-nghost+1][i]=H2[nghost+1][i];
                E1[nx1-nghost+1][i]=E2[nghost+1][i];
                T1[nx1-nghost+1][i]=T2[nghost+1][i];

            } break;

        case 'D':
            // Domain 2 is at the bottom of Domain 1
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx1-nghost; ++i)
            {
                rho1[i][0]=rho2[i][ny2-nghost-2];
                u1[i][0]=u2[i][ny2-nghost-2];
                v1[i][0]=v2[i][ny2-nghost-2];
                p1[i][0]=p2[i][ny2-nghost-2];
                H1[i][0]=H2[i][ny2-nghost-2];
                E1[i][0]=E2[i][ny2-nghost-2];
                T1[i][0]=T2[i][ny2-nghost-2];

                rho1[i][1]=rho2[i][ny2-nghost-1];
                u1[i][1]=u2[i][ny2-nghost-1];
                v1[i][1]=v2[i][ny2-nghost-1];
                p1[i][1]=p2[i][ny2-nghost-1];
                H1[i][1]=H2[i][ny2-nghost-1];
                E1[i][1]=E2[i][ny2-nghost-1];
                T1[i][1]=T2[i][ny2-nghost-1];
            } break;

        case 'U':
            // Domain 2 is at the top of Domain 1
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx1-nghost; ++i)
            {
                rho1[i][ny1-nghost]=rho2[i][nghost];
                u1[i][ny1-nghost]=u2[i][nghost];
                v1[i][ny1-nghost]=v2[i][nghost];
                p1[i][ny1-nghost]=p2[i][nghost];
                H1[i][ny1-nghost]=H2[i][nghost];
                E1[i][ny1-nghost]=E2[i][nghost];
                T1[i][ny1-nghost]=T2[i][nghost];

                rho1[i][ny1-nghost+1]=rho2[i][nghost+1];
                u1[i][ny1-nghost+1]=u2[i][nghost+1];
                v1[i][ny1-nghost+1]=v2[i][nghost+1];
                p1[i][ny1-nghost+1]=p2[i][nghost+1];
                H1[i][ny1-nghost+1]=H2[i][nghost+1];
                E1[i][ny1-nghost+1]=E2[i][nghost+1];
                T1[i][ny1-nghost+1]=T2[i][nghost+1];
            } break;

        default: printf("\nERROR: The specified location keyword is not recognized: '%c'!",location); exit(0); break;
    }
}


/* Slip-wall condition for a given face of a given fluid domain. At a slip wall, non-penetration of the flow is assumed by taking ghost cell normal velocity opposed to flow normal velocity. Other parameters are unchanged (zero gradient assumption). */
void slip_wall_condition(int nx,int ny,double **rho,double **p,double **u,double **v,double **H,double **E,double **T,char location,int nghost)
{   
    switch(location)
    {
        case 'L':
            // Wall condition on the left of domain
            #pragma omp parallel for schedule(static)
            for (int j = nghost; j < ny-nghost; ++j)
            {
                // First ghost-cell column ("0"-1)
                rho[1][j]=rho[2][j];
                p[1][j]=p[2][j];
                u[1][j]=-u[2][j];
                v[1][j]=v[2][j];
                H[1][j]=H[2][j];
                E[1][j]=E[2][j];
                T[1][j]=T[2][j];
                // Second ghost-cell column ("0"-2)
                rho[0][j]=rho[3][j];
                p[0][j]=p[3][j];
                u[0][j]=-u[3][j];
                v[0][j]=v[3][j];
                H[0][j]=H[3][j];
                E[0][j]=E[3][j];
                T[0][j]=T[3][j];
            } 
        break;

        case 'R':
            // Wall condition on the right...
            #pragma omp parallel for schedule(static)
            for (int j = nghost; j < ny-nghost; ++j)
            {
                // First ghost-cell column (N+1)
                rho[nx-2][j]=rho[nx-3][j];
                p[nx-2][j]=p[nx-3][j];
                u[nx-2][j]=-u[nx-3][j];
                v[nx-2][j]=v[nx-3][j];
                H[nx-2][j]=H[nx-3][j];
                E[nx-2][j]=E[nx-3][j];
                T[nx-2][j]=T[nx-3][j];
                // Second ghost-cell column (N+2)
                rho[nx-1][j]=rho[nx-4][j];
                p[nx-1][j]=p[nx-4][j];
                u[nx-1][j]=-u[nx-4][j];
                v[nx-1][j]=v[nx-4][j];
                H[nx-1][j]=H[nx-4][j];
                E[nx-1][j]=E[nx-4][j];
                T[nx-1][j]=T[nx-4][j];
            } 
        break;

        case 'D':
            // Wall condition at the bottom...
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx-nghost; ++i)
            {
                // First ghost-cell column ("0"-1)
                rho[i][1]=rho[i][2];
                p[i][1]=p[i][2];
                u[i][1]=u[i][2];
                v[i][1]=-v[i][2];
                H[i][1]=H[i][2];
                E[i][1]=E[i][2];
                T[i][1]=T[i][2];
                // Second ghost-cell column ("0"-2)
                rho[i][0]=rho[i][3];
                p[i][0]=p[i][3];
                u[i][0]=u[i][3];
                v[i][0]=-v[i][3];
                H[i][0]=H[i][3];
                E[i][1]=E[i][2];
                T[i][1]=T[i][2];
            } 
        break;

        case 'U':
            // Wall condition at the top...
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx-nghost; ++i)
            {
                // First ghost-cell column (N-1)
                rho[i][ny-2]=rho[i][ny-3];
                p[i][ny-2]=p[i][ny-3];
                u[i][ny-2]=u[i][ny-3];
                v[i][ny-2]=-v[i][ny-3];
                H[i][ny-2]=H[i][ny-3];
                E[i][ny-2]=E[i][ny-3];
                T[i][ny-2]=T[i][ny-3];
                // Second ghost-cell column (N-2)
                rho[i][ny-1]=rho[i][ny-4];
                p[i][ny-1]=p[i][ny-4];
                u[i][ny-1]=u[i][ny-4];
                v[i][ny-1]=-v[i][ny-4];
                H[i][ny-1]=H[i][ny-4];
                E[i][ny-1]=E[i][ny-4];
                T[i][ny-1]=T[i][ny-4];
            } 
        break;

        default: printf("\nERROR: The specified boundary location '%c' is not recognized!",location); exit(0); break;
    }
}


/* Zero-velocity wall condition. Is suitable only if the mesh is sufficiently fine at walls and if viscosity effects are taken into account. */
void noslip_wall_condition(int nx,int ny,double **rho,double **p,double **u,double **v,double **H,double **E,double **T,char location,int nghost)
{   
    switch(location)
    {
        case 'L':
            // Wall condition on the left of domain
            #pragma omp parallel for schedule(static)
            for (int j = nghost; j < ny-nghost; ++j)
            {
                // First ghost-cell column ("0"-1)
                rho[1][j]=rho[2][j];
                p[1][j]=p[2][j];
                u[1][j]=-u[2][j];
                v[1][j]=-v[2][j];
                H[1][j]=H[2][j];
                E[1][j]=E[2][j];
                T[1][j]=T[2][j];
                // Second ghost-cell column ("0"-2)
                rho[0][j]=rho[3][j];
                p[0][j]=p[3][j];
                u[0][j]=-u[3][j];
                v[0][j]=-v[3][j];
                H[0][j]=H[3][j];
                E[0][j]=E[3][j];
                T[0][j]=T[3][j];
            } 
        break;

        case 'R':
            // Wall condition on the right...
            #pragma omp parallel for schedule(static)
            for (int j = nghost; j < ny-nghost; ++j)
            {
                // First ghost-cell column (N+1)
                rho[nx-2][j]=rho[nx-3][j];
                p[nx-2][j]=p[nx-3][j];
                u[nx-2][j]=-u[nx-3][j];
                v[nx-2][j]=-v[nx-3][j];
                H[nx-2][j]=H[nx-3][j];
                E[nx-2][j]=E[nx-3][j];
                T[nx-2][j]=T[nx-3][j];
                // Second ghost-cell column (N+2)
                rho[nx-1][j]=rho[nx-4][j];
                p[nx-1][j]=p[nx-4][j];
                u[nx-1][j]=-u[nx-4][j];
                v[nx-1][j]=-v[nx-4][j];
                H[nx-1][j]=H[nx-4][j];
                E[nx-1][j]=E[nx-4][j];
                T[nx-1][j]=T[nx-4][j];
            } 
        break;

        case 'D':
            // Wall condition at the bottom...
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx-nghost; ++i)
            {
                // First ghost-cell column ("0"-1)
                rho[i][1]=rho[i][2];
                p[i][1]=p[i][2];
                u[i][1]=-u[i][2];
                v[i][1]=-v[i][2];
                H[i][1]=H[i][2];
                E[i][1]=E[i][2];
                T[i][1]=T[i][2];
                // Second ghost-cell column ("0"-2)
                rho[i][0]=rho[i][3];
                p[i][0]=p[i][3];
                u[i][0]=-u[i][3];
                v[i][0]=-v[i][3];
                H[i][0]=H[i][3];
                E[i][1]=E[i][2];
                T[i][1]=T[i][2];
            } 
        break;

        case 'U':
            // Wall condition at the top...
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx-nghost; ++i)
            {
                // First ghost-cell column (N-1)
                rho[i][ny-2]=rho[i][ny-3];
                p[i][ny-2]=p[i][ny-3];
                u[i][ny-2]=-u[i][ny-3];
                v[i][ny-2]=-v[i][ny-3];
                H[i][ny-2]=H[i][ny-3];
                E[i][ny-2]=E[i][ny-3];
                T[i][ny-2]=T[i][ny-3];
                // Second ghost-cell column (N-2)
                rho[i][ny-1]=rho[i][ny-4];
                p[i][ny-1]=p[i][ny-4];
                u[i][ny-1]=-u[i][ny-4];
                v[i][ny-1]=-v[i][ny-4];
                H[i][ny-1]=H[i][ny-4];
                E[i][ny-1]=E[i][ny-4];
                T[i][ny-1]=T[i][ny-4];
            } 
        break;

        default: printf("\nERROR: The specified boundary location '%c' is not recognized!",location); exit(0); break;
    }
}


/* Supersonic inlet condition, only available for left boundaries. The inlet being supersonic, all characteristics point inside the domain and therefore all must be specified. Here, incoming Mach, temperature and pressure are specified. Normal velocity component is supposed to be 0 (hence four conditions are specified). */
void supersonic_inlet_condition(int nx,int ny,double **rho,double **p,double **u,double **v,double **H,double **E,double **T,char location,int nghost,double gamma,double r,double Tref,double Mref,double Pref)
{
    // This condition will only be used on left boundaries
    if (location=='L')
    {
        #pragma omp parallel for schedule(static)
        for (int j = nghost; j < ny-nghost; ++j)
        {
            for (int i = 0 ; i < nghost ; i++)
            {
                rho[i][j]=Pref/(r*Tref);
                p[i][j]=Pref;
                u[i][j]=Mref*sqrt(gamma*r*Tref);
                v[i][j]=0.0;
                E[i][j]=p[i][j]/(gamma-1) + 0.5*rho[i][j]*pow(Mref*sqrt(gamma*p[i][j]/rho[i][j]),2);
                H[i][j]=(E[i][j]+p[i][j])/rho[i][j];
                T[i][j]=Tref;
            }
        }
    }
    else
    {
        fprintf(stderr,"\nERROR: Wrong direction specified as supersonic inlet condition!\n"); exit(0);
    }
}


/* Supersonic outlet condition only available at right and top boundaries. The outlet being supersonic, all characteristics point outside of the domain and so no additional condition needs to be specified. Here, zero gradient is assumed for all variables. */
void supersonic_outlet_condition(int nx,int ny,double **rho,double **p,double **u,double **v,double **H,double **E,double **T,char location,int nghost,double gamma,double r,double Tref,double Mref,double Pref)
{
    switch(location)
    {
        case 'R':
            #pragma omp parallel for schedule(static)
            for (int j = nghost; j < ny-nghost; ++j)
            {
                for (int i = 0 ; i < nghost ; i++)
                {
                    p[nx-1-i][j]=p[nx-nghost-1][j];
                    u[nx-1-i][j]=u[nx-nghost-1][j];
                    v[nx-1-i][j]=v[nx-nghost-1][j];
                    H[nx-1-i][j]=H[nx-nghost-1][j];
                    rho[nx-1-i][j]=rho[nx-nghost-1][j];
                    E[nx-1-i][j]=E[nx-nghost-1][j];
                    T[nx-1-i][j]=T[nx-nghost-1][j];
                }
            } break; 

        case 'U':
            #pragma omp parallel for schedule(static)
            for (int i = nghost; i < nx-nghost; ++i)
            {
                for (int j = 0 ; j < nghost ; j++)
                {
                    p[i][ny-1-j]=p[i][ny-nghost-1];
                    u[i][ny-1-j]=u[i][ny-nghost-1];
                    v[i][ny-1-j]=v[i][ny-nghost-1];
                    H[i][ny-1-j]=H[i][ny-nghost-1];
                    rho[i][ny-1-j]=rho[i][ny-nghost-1];
                    E[i][ny-1-j]=E[i][ny-nghost-1];
                    T[i][ny-1-j]=T[i][ny-nghost-1];
                }
            } break; 

        default: fprintf(stderr,"\nERROR: Wrong direction specified as supersonic outlet condition!\n"); exit(0); break;
    }
}


/* General function to update ghost cells between two time steps depending on the boundary condition that is specified there. If the user defines new boundary conditions, the latter should be mentionned also in this function so that the latter can be recognized. */
void update_ghost_cells(int Ndom,int *Nx,int *Ny,double ***x,double ***y,double ***xold,double ***yold,double ***rho,double ***u,double ***v,double***p,double ***H,double ***E,double ***T,char *b_loc,char ***b_type,double R_closed,int nghost,double gamma,double r,double Pref,double Tref,double Mref)
{
    #pragma omp parallel for
    for (int k = 0; k < Ndom; ++k)
    {
        #pragma omp parallel for num_threads(2)
        for (int jboundary = 0; jboundary < 4; ++jboundary)
        {
            // Slip-wall condition.
            if (strcmp(b_type[k][jboundary],"slip")==0)
            {
                slip_wall_condition(Nx[k],Ny[k],rho[k],p[k],u[k],v[k],H[k],E[k],T[k],b_loc[jboundary],nghost);
            }
            // No-slip wall condition.
            else if (strcmp(b_type[k][jboundary],"wall")==0)
            {
                noslip_wall_condition(Nx[k],Ny[k],rho[k],p[k],u[k],v[k],H[k],E[k],T[k],b_loc[jboundary],nghost);
            }
            // Connection between two adjacent fluid domains.
            else if (strncmp(b_type[k][jboundary],"con",3)==0)
            {
                int lengthString = strlen(b_type[k][jboundary]);
                char *lastChar = &b_type[k][jboundary][lengthString-1];
                int boundIndex = *lastChar - '0';
                connected_condition(Nx[k],Nx[boundIndex],Ny[k],Ny[boundIndex],rho[k],rho[boundIndex],p[k],p[boundIndex],u[k],u[boundIndex],v[k],v[boundIndex],H[k],H[boundIndex],E[k],E[boundIndex],T[k],T[boundIndex],b_loc[jboundary],nghost);
            }
            else if (strcmp(b_type[k][jboundary],"supI")==0)
            {
                supersonic_inlet_condition(Nx[k],Ny[k],rho[k],
                    p[k],u[k],v[k],H[k],E[k],
                    T[k],b_loc[jboundary],nghost,GAMMA,R,Tref,Mref,Pref);
            }
            else if (strcmp(b_type[k][jboundary],"supO")==0)
            {
                supersonic_outlet_condition(Nx[k],Ny[k],rho[k],
                    p[k],u[k],v[k],H[k],E[k],T[k],
                    b_loc[jboundary],nghost,GAMMA,R,Tref,Mref,Pref);
            }
            else
            {
                fprintf(stderr,"\nERROR: Boundary conditions specified were not recognized at domain %d: %s detected!\n",k,b_type[k][jboundary]);
            }
        }
    }
}

/* End of file */