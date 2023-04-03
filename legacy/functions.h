
/*  Collection of general functions used in the code. */

/* Mean value of any parameter stored in a matrix along a vertical or horizontal line.
    data: (nx,ny) sized matrix containing data to average
    nx,ny: Matrix dimensions
    nghost: Number of ghost cells at each boundary
    direction: V for vertical, H for horizontal
    index: Index at which to average data along the wanted direction
*/
double mean_along_line(double **data,int nx, int ny,int nghost,char direction,int index)
{
    double mean=0;

    switch(direction)
    {
        // Data should be averaged in the vertical direction
        case 'V': 
            for (int j = nghost; j < ny-nghost; ++j)
            {
                mean = mean + data[index][j]/(ny-2*nghost);
            }
        break;

        // Data should be averaged in the horizontal direction
        case 'H':
            for (int i = nghost; i < nx-nghost; ++i)
            {
                mean = mean + data[i][index]/(nx-2*nghost);
            }
        break;

        default: 
            printf("ERROR: Wrong specification in 'mean_along_line': %c!",direction); exit(0);
        break;
    }
    return mean;
}

/* Takes all the values in a rectangle around the valve, and returns the average value between all of them to approximate mean value around valve */
double mean_at_valve(int nx_start, int nx, int ny_start, int ny, double **data)
{
    double mean = 0.0;
    for (int i = nx_start; i < nx_start + nx; ++i)
    {
        for (int j = ny_start; j < ny_start + ny; ++j)
        {
            mean = mean + data[i][j]/ny/nx;
        }
    }
    return mean;
}

/* To compute the value of mass flow rate based on 3D CFD results using fitting functions and average values at the valve. */
double compute_mfr(double p_inf, double p_sup, double rho_sup, double y_tip, double l, double b0, double b1, double gamma, double r)
{
    double mfr_out;
    if (y_tip == 0.0)
    {
        mfr_out = 0.0;
    }
    else if (y_tip > 0.0)
    {
        // Reference area by Fukunari, pressure ratio and critical pressure ratio
        double a_ref = y_tip*(b1 + sqrt(l*l + 0.5*pow(b0-b1,2)));
        double pratio = fmax(0.0,fmin(p_inf/p_sup,1.0));
        double p_crit = pow(2.0/(gamma+1.0),gamma/(gamma-1.0));
        
        // OpenFOAM fitting of discharge coefficient (Cd) is specified here!
        // The fit function (equation 2.3 in Florian (2017)) based on cfd simulations to get a discharge coefficient.
        double cd_3d = fmax(0.0,fmin(1.0,0.9193*pow(1.0+0.5*y_tip*1000.0,-0.596)));
        // printf("CD3D: %f\n", cd_3d);

        if (pratio <= 0.0 || pratio > 1.0)
        {
            mfr_out = 0.0;
        }
        else if (pratio > 0.0 && pratio <= p_crit)
        {
            mfr_out = cd_3d*a_ref*gamma*p_sup/sqrt(gamma*p_sup/rho_sup)*pow(2.0/(gamma+1.0),0.5*(gamma+1.0)/(gamma-1.0));
        }
        else if (pratio > p_crit && pratio <= 1.0)
        {
            mfr_out = cd_3d*a_ref*rho_sup*pow(pratio,1.0/gamma)*sqrt(2.0*gamma*p_sup/rho_sup/(gamma-1.0)*(1.0-pow(pratio,(gamma-1.0)/gamma)));
        }
        
        // printf("MFR IS: %f\n",*mfr_out);
        if (isnan(p_sup) || isnan(rho_sup))
        {
            fprintf(stderr, "\nERROR: RHOSUP OR PSUP IS NaN: %f %f!\n",p_sup,rho_sup);
            exit(0);
        }
        if (isnan(mfr_out))
        {
            fprintf(stderr, "\nERROR: MFR COMPUTED IS NaN: PSUP AND PINF ARE %f %f\n",p_sup,p_inf);
        exit(0);
        }
    }
    else if (y_tip < 0.0)
    {
        fprintf(stderr, "\nERROR: Valve displacement is negative: MFR cannot be computed! Y_TIP = %g!\n",y_tip);
        exit(0);
    }

    return mfr_out;
}

/* Compute the total source term based on mass flow rate and fluid variables around the reed valve.
Sets the source[4] by reference, where the index corresponds to the euler equation it will be coupled to: [density term, x-axis velocity, y/r-axis velocity, and specific internal energy]

*/ 
void compute_cell_source(int k_dom, int k_inf, int k_sup, double r_low, double r_up, double mfr_valve, double y_tip, double theta_tip, double dy, int n_per_stage, double gamma,double radius, double l_hole, double hole_factor, double dx, double b_hole, int mfr_n_cell,double rho_mean_sup, double p_mean_sup,double p_mean_inf, double rho_sup, double u_sup, double v_sup, double p_sup, double b0, double b1, double l, double source[4])
{   
    source[0] = 0.0;
    source[1] = 0.0;
    source[2] = 0.0;
    source[3] = 0.0;

    if (y_tip > 0.0 && mfr_valve > 0.0 && p_mean_inf < p_mean_sup) // threshold value of tip displacement
    {
        // Equivalent surface for one interface cell (applied to velocity)
        // double A_eq = l_hole/2.0*(2.0*M_PI*radius)/mfr_n_cell;
        // double A_eq = n_per_stage*l_hole/2.0*b_hole/mfr_n_cell;
        // Same but applied to pressure
        // double A_hole = n_per_stage*l_hole*b_hole/mfr_n_cell;

        // The fit function (equation 2.3 in Florian (2017)) based on cfd simulations to get a discharge coefficient.
        // This function call is only referred to by the area, which itself is not used- remove?
        double cd_3d = fmax(0.0,fmin(1.0,0.9193*pow(1.0+0.5*y_tip*1000.0,-0.596)));
        // double A_eq = cd_3d*n_per_stage*y_tip*(b1 + sqrt(l*l + 0.5*pow(b0-b1,2)));

        // This equivalent area is never used!
        double A_eq = cd_3d*n_per_stage*l_hole*b_hole/mfr_n_cell;
        // double A_eq = n_per_stage*l_hole*b_hole*sin(theta_tip)/mfr_n_cell;

        double r_cell, v_cell;
        if (k_dom == k_inf)
        {
            r_cell = r_low;
            // v_cell = l_hole*hole_factor*(pow(r_cell,2)-pow(r_cell-dy,2));
            v_cell = dx*M_PI*(pow(r_cell,2)-pow(r_cell-dy,2));
        }
        else if (k_dom == k_sup)
        {
            r_cell = r_up;
            // v_cell = l_hole*hole_factor*(pow(r_cell+dy,2)-pow(r_cell,2));
            v_cell = dx*M_PI*(pow(r_cell+dy,2)-pow(r_cell,2));
        }

        // Total mass flow rate per cell
        // double mfr_eq = (n_per_stage*mfr_valve);
        double mfr_eq = (n_per_stage*mfr_valve)/mfr_n_cell;

        // Normal velocity assuming that total mass flow rate is eq. to rho*A_eq*V_eq
        double p_ratio = fmax(0.0,fmin(1.0,p_mean_inf/p_mean_sup));
        double p_crit = pow(2.0/(gamma+1.0),gamma/(gamma-1.0));
        double rho_eq = 0.0, p_eq = 0.0;
        if (p_ratio > 0.0 && p_ratio <= p_crit) // choked flow
        {
            rho_eq = rho_mean_sup*pow(p_crit,1.0/gamma);
            p_eq = p_crit*p_mean_sup;
        }
        else if (p_ratio > p_crit && p_ratio <= 1.0) // non-choked flow
        {
            rho_eq = rho_mean_sup*pow(p_ratio,1.0/gamma);
            p_eq = p_mean_inf;
        }
        else
        {
            fprintf(stderr, "\nERROR: UNEXPECTED SOURCE TERM PRESSURE RATIO: %f %f %f",p_mean_inf,p_mean_sup,p_ratio);
            exit(0);
        }

        rho_eq=rho_mean_sup;
        p_eq=p_mean_sup;
        // double V_eq = mfr_eq/rho_eq/A_eq;
        // double V_eq = mfr_eq/rho_eq/(2.0*M_PI*r_cell*l_hole*hole_factor);
        double A_cyl = (2.0*M_PI*r_cell*dx);
        double V_eq = mfr_eq/rho_eq/A_cyl;
        double source_u = (mfr_eq*u_sup)/v_cell;
        double source_v = (mfr_eq*v_sup)/v_cell;
        // double source_v = (mfr_eq*V_eq )/v_cell;

        // printf("V_EQ: %f\n", V_eq);
        
        // Source for lower domain
        if (k_dom == k_inf)
        {
            source[0] = mfr_eq/v_cell;
            source[1] = source_u;// + cos(theta_tip)*(p_eq)*A_eq/v_cell;
            source[2] = source_v;// + (p_eq)*A_eq/v_cell;
            source[3] = (gamma/(gamma-1.0)*p_eq/rho_eq + 0.5*V_eq*sqrt(u_sup*u_sup+v_sup*v_sup))*mfr_eq/v_cell;
            if (isnan(source[0]) || isnan(source[1]) || isnan(source[2]) || isnan(source[3]))
            {
                printf("%f\n",source[0]);
                printf("%f\n",source[1]);
                printf("%f\n",source[2]);
                printf("%f\n",source[3]);
            }
        }
        // Source for upper domain
        else if (k_dom == k_sup)
        {
            source[0] = -mfr_eq/v_cell;
            source[1] = -source_u;
            source[2] = -source_v;// - (p_eq)*A_eq/v_cell;
            source[3] = -(gamma/(gamma-1.0)*p_eq/rho_eq + 0.5*V_eq*sqrt(u_sup*u_sup+v_sup*v_sup))*mfr_eq/v_cell;
            if (isnan(source[0]) || isnan(source[1]) || isnan(source[2]) || isnan(source[3]))
            {
                printf("%f\n",mfr_eq);
                printf("%f\n",u_sup);
                // printf("%f\n",V_eq);
                printf("%f\n",p_sup/rho_sup);
            }
        }

    }    
}

/* Allocates and initializes an array of size N or (N+1) depending on gridOrCell argument.
    N: Length of array (typically, number of cells)
    gridOrCell: 1 to add an aditionnal element (N cells will require N+1 nodes)
*/
double * init_array(int N,int gridOrCell)
{
    double * data = malloc((N+gridOrCell)*sizeof(double));
    if(data == NULL)
    {
        printf("\nERROR: ERROR OF ALLOCATION IN 'init_matrix'!");
    }

    for (int i = 0; i < N+gridOrCell; ++i)
    {
        // All coefficients are initialized to 0.0
        data[i]=0.0;
    }
    return data;
}

/* Allocates and initializes a matrix of size NxM or (N+1)x(M+1) depending on gridOrCell argument 
    Nx,Ny: Size of matrix
    gridOrCell: Same as before
*/
double ** init_matrix(int Nx, int Ny,int gridOrCell)
{
    double ** data = malloc((Nx+gridOrCell)*(Ny+gridOrCell)*sizeof(double *));
    if(data == NULL)
    {
        printf("\nERROR: ERROR OF ALLOCATION IN 'init_matrix'!");
    }

    // For each line...
    for (int i = 0; i < Nx+gridOrCell; ++i)
    {
        // Each line has the size of its respective number of columns Ny+gridOrCell
        data[i] = malloc((Ny+gridOrCell)*sizeof(double));
        if ( data[i]==NULL )
        {
            printf("\nERROR: Out of memory in matrix-pool allocation (3)!");
        }
        // For each column...
        for (int j = 0; j < Ny+gridOrCell; ++j)
        {
            // All coefficients are initialized to 0.0
            data[i][j]=0.0;
        }
    }
    return data;
}

/* Generate an int matrix */
int ** init_int_matrix(int Nx, int Ny, int gridOrCell)
{
    int ** data = malloc((Nx+gridOrCell)*(Ny+gridOrCell)*sizeof(int *));
    if(data == NULL)
    {
        printf("\nERROR: ERROR OF ALLOCATION IN 'init_matrix'!");
    }

    // For each line...
    for (int i = 0; i < Nx+gridOrCell; ++i)
    {
        // Each line has the size of its respective number of columns Ny+gridOrCell
        data[i] = malloc((Ny+gridOrCell)*sizeof(int));
        if ( data[i]==NULL )
        {
            printf("\nERROR: Out of memory in matrix-pool allocation (3)!");
        }
        // For each column...
        for (int j = 0; j < Ny+gridOrCell; ++j)
        {
            // All coefficients are initialized to 0.0
            data[i][j]=0;
        }
    }
    return data;
}

/* Initialize the arrays representing variables for each cell in all domains, all at once.
    Ndom,Nx[],Ny[]: Same as before
    rho,u,v,p,E,T,H: Variables to initialize*/
void init_domain(int Ndom,int Nx[Ndom],int Ny[Ndom],double ****rho,double ****u,double ****v,double ****p,double ****E,double ****T,double ****H)
{
    // Fluid variables
    *rho = init_variable(Ndom,Nx,Ny,0);
    *u = init_variable(Ndom,Nx,Ny,0);
    *v = init_variable(Ndom,Nx,Ny,0);
    *p = init_variable(Ndom,Nx,Ny,0);
    *E = init_variable(Ndom,Nx,Ny,0);
    *T = init_variable(Ndom,Nx,Ny,0);
    *H = init_variable(Ndom,Nx,Ny,0);
}

double ** init_interface(int n_valve, int *n_cell_interface)
{
    int total_length = 0;
    for (int i = 0; i < n_valve; ++i)
    {
        total_length += n_cell_interface[i];
        // printf("%i %i\n", n_cell_interface[i],total_length);
    }

    double **array = malloc(total_length*sizeof(double*));
    for (int i = 0; i < n_valve; ++i)
    {
        array[i] = malloc(n_cell_interface[i]*sizeof(double));
        for (int j = 0; j < n_cell_interface[i]; ++j)
        {
            array[i][j] = 0.0;
        }
    }

    return array;
}

/* Looks for highest velocity in computational domain and computes CFL accordingly.
    u,v,T: velocity and temperature data used to estimate Mach number and CFL
    Ndom,Nx[],Ny[]: Same as before
    dx,dy: Grid sizes along x and y
    dt: Time step
*/
double get_cfl(double ***x,double ***y,double ***u,double ***v,double ***T, double r,double gamma,int Ndom,int Nx[Ndom],int Ny[Ndom],double dt)
{
    double CFL=0;

    #pragma omp parallel for schedule(static) reduction(max:CFL)
    for (int k = 0; k < Ndom; ++k)
    {
        double CFL_dom=0.0;

        #pragma omp parallel for reduction(max:CFL_dom)
        for (int i = 0; i < Nx[k]; ++i)
        {
            double umax,vmax,dx,dy;

            for (int j = 0; j < Ny[k]; ++j)
            {
                dx = x[k][i+1][j] - x[k][i][j];
                dy = y[k][i][j+1] - y[k][i][j];
                umax = fabs(u[k][i][j])+sqrt(gamma*r*T[k][i][j]);
                vmax = fabs(v[k][i][j])+sqrt(gamma*r*T[k][i][j]);
                // a = sqrt(gamma*r*T[k][i][j]);
                // CFL_dom = fmax(CFL_dom,dt*(umax/dx + vmax/dy + a*sqrt(1.0/pow(dx,2)+1.0/pow(dy,2))));
                // CFL_dom = fmax(CFL_dom,dt*sqrt(umax*umax/dx/dx+vmax*vmax/dy/dy));
                CFL_dom = fmax(CFL_dom,dt*fmax(umax/dx,vmax/dy));
            }
        }

        CFL = CFL_dom;
    }

    return CFL;
}

/* Theoretical value of supersonic stagnation pressure */
double sup_stag_p(double m0,double p0,double gamma)
{
    double var = p0*pow((gamma+1.0)*m0*m0/2.0,gamma/(gamma-1.0))*pow((gamma+1.0)/(2.0*gamma*m0*m0-gamma+1.0),1.0/(gamma-1.0));
    return var;
}

/* Theoretical value of supersonic stagnation temperature */
double sup_stag_t(double m0,double t0,double gamma)
{
    // double m2 = sqrt((1.0+(gamma-1.0)/2.0*m0*m0)/(gamma*m0*m0-(gamma-1.0)/2.0));
    // double t2 = t0*(2.0*gamma*m0*m0-(gamma-1.0))*(2.0+(gamma-1.0)*m0*m0)/pow((gamma+1.0)*m0,2);
    double t2 = t0*(1.0+(gamma-1.0)/2.0*m0*m0);
    // double ti2 = t2*(1.0+(gamma-1.0)*m2*m2/2.0);
    // return ti2;
    return t2;
}

/* Import a matrix from an existing .DAT file */
void import_matrix(int k, char *var, int nx, int ny, int nghost, double **matrix)
{
    char index[5],filename[30];
    sprintf(index,"%i",k);
    strcpy(filename,PRERUN_FOLDER);
    strcat(filename,index);
    strcat(filename,"/");
    strcat(filename,var);
    strcat(filename,EXP_EXTENSION);
    // printf("%s\n",filename);

    FILE *fp = fopen(filename,"r");

    for (int i = nghost; i < nx-nghost; ++i)
    {
        for (int j = nghost; j < ny-nghost; ++j)
        {
            if (!fscanf(fp,"%lf",&matrix[i][j]))
            break;
        }
    }
}

/* Import the average of a given field for a given subdomain of any pre-run simulation */
double import_average(int k, char *var, int nx, int ny, int nghost)
{
    char index[5],filename[30];
    sprintf(index,"%i",k);
    strcpy(filename,PRERUN_FOLDER);
    strcat(filename,index);
    strcat(filename,"/");
    strcat(filename,var);
    strcat(filename,EXP_EXTENSION);
    // printf("%s\n",filename);

    FILE *fp = fopen(filename,"r");
    double number = 0.0;
    double average = 0.0;

    for (int i = nghost; i < nx-nghost; ++i)
    {
        for (int j = nghost; j < ny-nghost; ++j)
        {
            if (!fscanf(fp,"%lf",&number))
            break;
            average += number/(nx-2*nghost)/(ny-2*nghost);
        }
    }
    // printf("Average is : %f\n",average);

    return average;
}

// double domain_average(int nx, int ny, double **x, double **r, double **data, int nghost)
double domain_average(int nx, int ny, double **x, double **y, double **data, int nghost)
{
    // r is the y coordinate (cylindrical assumption)
    double value=0.0;
    double volume = (M_PI)*(pow(y[nghost][ny-nghost],2)-pow(y[nghost][nghost],2))*(x[nx-nghost][nghost]-x[nghost][nghost]);
    // double r_min = r[nghost][nghost];
    // double r_sup = r[nghost][ny-nghost];
    // double l = x[nx-nghost][nghost] - x[nghost][nghost];
    // p[dom_low][NGHOST][i]*(pow(r[dom_low][NGHOST][i+1],2)-pow(r[dom_low][NGHOST][i],2))/pow(R0,2);
    for (int i = nghost; i < nx-nghost; ++i)
    {
        for (int j = nghost; j < ny-nghost; ++j)
        {
            // value += data[i][j]/(nx-2*nghost)/(ny-2*nghost);
            value += data[i][j]*M_PI*fabs(pow(y[i][j+1],2)-pow(y[i][j],2))*(x[i+1][j]-x[i][j])/volume;
        }
    }

    return value;
}

// Calculate the new V(normal) and U(parallel) velocities with regard to a wall inclined from theta from the horizontal direction
void theta_frame_velocity(double u_x,double v_y, double theta)
{
    double u_theta = u_x*cos(theta)+v_y*sin(theta);
    double v_theta = -u_x*sin(theta)+v_y*sin(theta);
    u_x = u_theta;
    v_y = v_theta;
}

/* Calculate mass flow rate at a given face of any subdomain (only valid in the horizontal direction) */
double mfr_face(int nx, int ny, double **x, double **y, double **rho, double **u, int nx0, int nghost)
{
    double mfr = 0.0;
    for (int j = nghost; j < ny-nghost; ++j)
    {
        mfr += rho[nx0][j]*u[nx0][j]*M_PI*fabs(pow(y[nx0][j+1],2)-pow(y[nx0][j],2));
    }
    return mfr;
}

/* Face average for any subdomain (only for faces whose surface normal vector is horizontal) */ 
double average_face(int nx, int ny, double **x, double **y, double **p, int nx0, int nghost)
{
    double ave = 0.0;
    for (int j = nghost; j < ny-nghost; ++j)
    {
        ave += p[nx0][j]*M_PI*fabs(pow(y[nx0][j+1],2)-pow(y[nx0][j],2));
    }
    ave = ave/(M_PI)/(pow(y[nx0][ny-nghost],2)-pow(y[nx0][nghost],2));

    return ave;
}


// END OF FILE