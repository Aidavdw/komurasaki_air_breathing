/*  TITLE: 2D C# MICROWAVE BEAMED ENERGY ROCKET MODEL
    AUTHOR: Florian NGUYEN - The University of Tokyo, 2016-2017

    DESCRIPTION: 
    A 2D compressible, inviscid solver modelling a microwave rocket tube, several stages of reed valves, and a plenum region compressing inlet flow.
    Euler equations are solved by means of AUSM-DV scheme and piecewise parabolic MUSCL interpolation with Van Albada limiter (can be modified). 
    The reed valve is modelled with a Finite Element model as a double tapered, cross-section varying beam, based on which stiffness and mass matrices are computed. Structural damping is assumed to be of Rayleigh type (a*M + b*K), where a and b are fitted with experimental data (free vibration of the valve by Fukunari).
*/

#define getName(var)  #var
#define M_PI (4.0*atan(1.0))
#define _POSIX_C_SOURCE 199309L

#include <stdio.h>      // Basic C# header
#include <stdlib.h>     // Basic C# header
#include <math.h>       // Mathematical expressions
#include <string.h>     // To handle strings and char arrays
#include <time.h>       // For CPU time calculation
#include <omp.h>        // For multi-processing with Open MP
#include <sys/stat.h>   // To create directories

#include "legacy/colors.h"      // Colors for use in terminal
#include "legacy/parameters.h"  // List of parameters specified by user
#include "legacy/functions.h"   // General functions
#include "initial.h"     // Initial conditions suitable for each case
#include "legacy/muscl.h"       // MUSCL interpolation
#include "legacy/ausm.h"        // AUSM scheme
#include "microwave_orig.h"   // Microwave heating and MSD theory
#include "legacy/export.h"      // Export functions
#include "legacy/bound_cond.h"  // Boundary conditions
#include "legacy/fem_model.h"   // FEM model of the reed valve

#include "sim_case.h"


/* MAIN ROUTINE */
int main()
{
    /* STARTING PROGRAM, GETTING TIME... */
    setbuf(stdout, NULL);
    printf("\nProgram starting...\n");
    struct timespec start_time, end_time,start_time_loop,end_time_loop;
    clock_gettime( CLOCK_REALTIME, &start_time);

    // INITIALIZE CASE
    double L_T;         // Length of the reed valve; fixed + free part
	double M0;
	double * XSTART;
	double * YSTART;
	double * XLENGTH;
	double * YLENGTH;
	double * GRID_RATIO_X;
	double * GRID_RATIO_Y;
    double *X_V_START; // Place where each individual valve FEM cell starts (bottom left)
    char *B_LOC, ***B_TYPE;
    int* NXtot;         // The total amount of cells in the X-direction. Is an array, index is domain.
    int* NYtot;         // The total amount of cells in the Y-direction. Is an array, index is domain.
    int NDOMAIN = 0;    // The amount of domains in the simulation case.
    int N_VALVE = 0;    // The amount of valves used for the simulation
    int SOLID_ON = 0;   // I think this is a flaf on whether or not the FEM modeling should be used for calculating the deflection of the reed valves.
    int dom_low;        // I think this is the id of the domain that is on the bottom; the one that actually has the tube.
    int dom_up;         // I think this is the id of the domain that is on top; the one that feeds the valves.
    int PLENUM_ON = 0;  // I think this is a flag that is set on whether or not to use the plenum
    init_case(SIM_CASE,&L_T,&M0,&dom_low,&dom_up,&XSTART,&YSTART,&XLENGTH,&YLENGTH,&GRID_RATIO_X,&GRID_RATIO_Y,&X_V_START,&B_LOC,&B_TYPE,&NXtot,&NYtot,&NDOMAIN,&N_VALVE,&SOLID_ON,&PLENUM_ON);

    /* TIME STEP, EXPORT SPAN AND BASIC GRID SIZE */
    int Ntstep = (int)(1+TSIM/DT);                           // Number of time steps
    int N_CFL = (int)fmax(1,(int)(Ntstep)/N_STEP_CFL);       // Iterations between two displays of CFL
    int N_EXPORT = (int)fmax(1,(int)(Ntstep)/N_STEP_EXP);    // Iterations between two exports of data


    /* CREATE OUTPUT FOLDER */
    make_dir(OUT_FOLDERNAME);


    /* TOTAL NUMBER OF CELLS IN CURRENT SIMULATION, ORDER OF DX/DY */
    int Ncell = 0;
    for (int i = 0; i < NDOMAIN; ++i)
    {
        Ncell += (NXtot[i]-2*NGHOST)*(NYtot[i]-2*NGHOST);
    }
	printf("Total number of cells (without ghost cells) is: " BLU "%d" RESET " cells.\n",Ncell);
    printf("Reference cell size based on domain 0 is: DX="BLU "%f" RESET " m and DY=" BLU "%f" RESET " m.\n",XLENGTH[0]/(NXtot[0]-2*NGHOST),YLENGTH[0]/(NYtot[0]-2*NGHOST));
    

    /* CHECKING DOMAIN BOUNDARIES */
    printf("\nVerification of boundaries for each domain:\n");
    for (int i = 0; i < NDOMAIN; ++i)
    {
        printf("DOMAIN %d (%d x %d cells) boundaries: %s - %s - %s - %s.\n",i,NXtot[i]-2*NGHOST,NYtot[i]-2*NGHOST,B_TYPE[i][0],B_TYPE[i][1],B_TYPE[i][2],B_TYPE[i][3]);
    }


    /* DISPLAY INFORMATION ON TERMINAL */
    printf("\nCFL will be displayed every " BLU "%d" RESET " iterations (%.2f ms).\n",N_CFL,N_CFL*DT*1000);
    printf("Results will be exported every " BLU "%d " RESET "iterations (%.2f ms).\n",N_EXPORT,N_EXPORT*DT*1000);
    printf("Total number of iterations is: " BLU "%d " RESET "iterations.\n",Ntstep-1);
    printf("(Pref,Tref,Mref) = " BLU "(%f, %f, %f)" RESET ".\n", P0,T0,M0);


    /* VARIABLE DEFINITION AND ALLOCATION FOR FLUID MODEL */
    double ***x;        // All the bottom-left x-positions for all cells
	double  ***y;       // All the bottom-left y-positions for all cells
	double  ***xc;      // the x coordinate of the centre for all cells
	double  ***yc;      // the x coordinate of the centre for all cells
	double  ***xold;
	double  ***yold;
    double ***p;        // the pressure for all cells 
	double  ***rho;     // the pressure for all cells 
	double  ***u;       // the pressure for all cells 
	double  ***v;       // the pressure for all cells 
	double  ***E;       // the internal energy for all cells
	double  ***T;       // the temperature for all cells
	double  ***H;       // the enthalpy for all cells 
    double ***pt;		// field variable that will be set as the main one for the next time step; the iteration buffer.
	double  ***rhot;	// field variable that will be set as the main one for the next time step; the iteration buffer.
	double  ***ut;		// field variable that will be set as the main one for the next time step; the iteration buffer.
	double  ***vt;		// field variable that will be set as the main one for the next time step; the iteration buffer.
	double  ***Et;		// field variable that will be set as the main one for the next time step; the iteration buffer.
	double  ***Tt;		// field variable that will be set as the main one for the next time step; the iteration buffer.
	double  ***Ht;		// field variable that will be set as the main one for the next time step; the iteration buffer.
    double*** pRK;      // runge kutta field buffer.
    double*** rhoRK;	// runge kutta field buffer.
    double ***uRK;		// runge kutta field buffer.
	double  *** vRK;	// runge kutta field buffer.
	double  *** ERK;	// runge kutta field buffer.
	double  *** TRK;	// runge kutta field buffer.
	double  *** HRK;	// runge kutta field buffer.
	double  *** sonic_x;	// mask corresponding to every cell. If > 1, MUSCL has identified that the flow goes supersonic in the x-direction. if 0, it's subsonic.
	double  *** sonic_y;	// mask corresponding to every cell. If > 1, MUSCL has identified that the flow goes supersonic in the y-direction. if 0, it's subsonic.
    init_domain(NDOMAIN,NXtot,NYtot,&rho,&u,&v,&p,&E,&T,&H);
    init_domain(NDOMAIN,NXtot,NYtot,&rhot,&ut,&vt,&pt,&Et,&Tt,&Ht);
    init_domain(NDOMAIN,NXtot,NYtot,&rhoRK,&uRK,&vRK,&pRK,&ERK,&TRK,&HRK);
    x = init_variable(NDOMAIN,NXtot,NYtot,1);
    y = init_variable(NDOMAIN,NXtot,NYtot,1);
    xold = init_variable(NDOMAIN,NXtot,NYtot,1);
    yold = init_variable(NDOMAIN,NXtot,NYtot,1);
    xc = init_variable(NDOMAIN,NXtot,NYtot,0);
    yc = init_variable(NDOMAIN,NXtot,NYtot,0);
    sonic_x = init_variable(NDOMAIN,NXtot,NYtot,0);
    sonic_y = init_variable(NDOMAIN,NXtot,NYtot,0);

    if (p==NULL || u==NULL || v==NULL || rho==NULL || T==NULL || H==NULL || pt==NULL || ut==NULL || vt==NULL || rhot==NULL || Tt==NULL || Ht==NULL || pRK==NULL || uRK==NULL || vRK==NULL || rhoRK==NULL || TRK==NULL || HRK==NULL || x==NULL || y==NULL || xc==NULL || yc==NULL)
    {
        fprintf(stderr, "\nERROR: out of memory during allocation! \n");
        exit(0);
    }
    printf("\nAll variables generated successfully...\n");

    /* MESHING ALL DOMAINS (NEW AND OLD ALIKE) */
    compute_mesh(NDOMAIN,NXtot,NYtot,x,y,xc,yc,XSTART,YSTART,XLENGTH,YLENGTH,GRID_RATIO_X,GRID_RATIO_Y,NGHOST);

    /* EXPORT PARAMETERS REQUIRED FOR POST-PROCESSING AND SOLUTION RECONSTRUCTION */
    export_parameters(NDOMAIN,TSIM,DT,XSTART,YSTART,XLENGTH,YLENGTH,N_VALVE,N_FEM,NXtot,NYtot,NGHOST,N_EXPORT,W_FORMAT,PAR_FILENAME);


    /* INITIAL CONDITIONS BASED ON MICROWAVE DETONATION THEORY */
    double M_MSD, P1, U1, RHO1, M1, P2, RHO2, L_EXP;
    double p_prerun, T_prerun;
    int N_ITER_MSD;
    // If prerun case is involved, replace the reference P/T couple for post-MSD calculation
    if (strcmp(SIM_CASE,"sup_plen_rocket")==0)
    {
        // rho_prerun = import_average(dom_up-3,"rho",Nxtot[dom_up-3],NYtot[dom_up-3],NGHOST);
        p_prerun = import_average(dom_up-3,"p",NXtot[dom_up-3],NYtot[dom_up-3],NGHOST);
        T_prerun = import_average(dom_up-3,"T",NXtot[dom_up-3],NYtot[dom_up-3],NGHOST);
        N_ITER_MSD = solve_MSD(&M_MSD, &P1, &U1, &RHO1, &M1, &P2, &RHO2, &L_EXP,T_prerun,p_prerun,ETA,S0,R,GAMMA,L_TUBE,R0);
        printf("\nPRERUN CONDITIONS ARE: \nP=%f, \nT=%f. \n",p_prerun,T_prerun);
    }
    else
    {
        N_ITER_MSD = solve_MSD(&M_MSD, &P1, &U1, &RHO1, &M1, &P2, &RHO2, &L_EXP,T0,P0,ETA,S0,R,GAMMA,L_TUBE,R0);
    }
    printf("MSD conditions in the rocket were estimated after %d iterations.\n",N_ITER_MSD);


    /* INITIAL CONDITIONS ON DOMAINS */
    apply_initial_conditions(SIM_CASE,NDOMAIN,NXtot,NYtot,x,rho,u,v,p,E,T,H,NGHOST,M_MSD,P1,U1,RHO1,M1,P2,RHO2,L_EXP,L_TUBE,T0,P0,M0,R,GAMMA);
    printf("\nInitial conditions applied...\n");

    /* EXPORT INITIAL SOLUTION AND CREATE OUTPUT FOLDER, DEFINE WALL PRESSURE AND INTAKE MASS FLOW RATE */
    export_fluid_data(NDOMAIN,NXtot,NYtot,x,y,xc,yc,rho,u,v,p,E,T,H,0.0,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT);
    double p_wall = p[dom_low][NGHOST][NGHOST];
    double mfr_intake = 0.0, p_tube_mean = 0.0, rho_tube_mean = 0.0;
    double p_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NGHOST);
    double rho_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],NGHOST);
    double p_tube = domain_average(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],p[dom_low],NGHOST);
    double rho_tube = domain_average(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],rho[dom_low],NGHOST);
    double mfr_plenum = mfr_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],u[dom_up],NGHOST,NGHOST);
    double p_drag = average_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NXtot[dom_up]-NGHOST-1,NGHOST);

    append_data_to_file(p_wall,0.0,P_MEAN_WALL_FILENAME,EXP_EXTENSION);
    append_data_to_file(p_tube,0.0,P_TUBE_FILENAME,EXP_EXTENSION);
    append_data_to_file(rho_tube,0.0,RHO_TUBE_FILENAME,EXP_EXTENSION);
    append_data_to_file(mfr_intake,0.0,INTAKE_MFR_FILENAME,EXP_EXTENSION);
    append_data_to_file(0.0,0.0,PFR_FILENAME,EXP_EXTENSION);
    append_data_to_file(0.0,0.0,MFR_TOT_FILENAME,EXP_EXTENSION);

    if (PLENUM_ON)
    {
        append_data_to_file(p_drag,0.0,PLENUM_DRAG_PRESSURE_FILENAME,EXP_EXTENSION);
        append_data_to_file(p_plenum,0.0,PLEN_P_FILENAME,EXP_EXTENSION);
        append_data_to_file(rho_plenum,0.0,PLEN_RHO_FILENAME,EXP_EXTENSION);
        append_data_to_file(mfr_plenum,0.0,PLENUM_MFR_FILENAME,EXP_EXTENSION);
    }

    printf("\nInitial fluid conditions exported...\n");

    /* DEFINITION AND ALLOCATION OF FEM VARIABLES */
	int N_DOF;              // Total amount of degrees of freedom for FEM values
	int  N_NODE;            // The amount of nodes that are in the FEM simulation. NFEM + 1
	int  N_ACTIVE;          // The amount of nodes that is allowed to deform. FemDeformation::freeNodes * DOF_PER_NODE
	int  N_INACTIVE;        // The amount of nodes that is forced to stay fixed in place. FemDeformation::fixedNodes * DOF_PER_NODE
	int  *act_DOF;          // A vector that just increases by one every step, going from (n_node * N_DOF) -> (n_total * DOF + DOF-1).
	double **x_FEM;         // X locations of the inidividual FEM cells, relative to the left-most point.
	double **y_FEM;			// A redundant extra copy of the locations of the individual fem cells. Set from u2. 
	double *b;              // The (average) width of the cell. Literally b in Florian (2017).
	double *h;
	double *A;
	double *I;
    double **K;             // Stiffness matrix. Size of dof*dof
	double **C;             // Structural damping matrix.
	double **M;             // global mass matrix for 1 valve FEM beam
	double **L_K;           // cholesky decomp. of the stiffness matrix.
	double **LT_K;          // transpose of the cholesky decomp. of the stiffness matrix. only used to calculate L_K
	double **L_R1;              // Cholesky decomp. of newmark matrix 1, used in newmark solving.
	double **LT_R1;             // Transpose of the Cholesky decomposition of newmark matrix one. only used to calculate L_R1
	double **R1;                // Newmark matrix 1. This one is only used for cholesky decomposition. L_R1 is the one that's actually used for solving.
	double **R2;                // Newmark matrix 2. Is used for newmark solving.
	double **R3;                // Newmark matrix 3
    double **U0_DOF;			// Some deflection of the reed valve? Input in update_valve() but is never actually referenced.
	double **U1_DOF;			// Deflection of the reed valve node in the previous time step. in update valve, is used if u2[i] > 0.012 to set the physical deflection. After initialisation, this is set to U2_DOF
	double **U2_DOF;			// Deflection of the reed valve node in the current time step. in update valve, is used if u2[i] < 0.012 to set the physical deflection. After initialisation, this is set to U2_DOF
	double **U2_DOF_K;			// Some deflection of the reed valve?
	double **dp_interface;		// Doesn't appear to be used.
	double **p_FEM;				// The difference in pressure between the FEM node's set sample cell and the ambient conditions in Pa. 
	double**F_DOF;				// Holds all the force components, to be used in solving the FEM system of equations.

    /* MASS-FLOW RATE RELATED PARAMETERS DUE TO VALVE */
    // To locate cells with source term
    int *mfr_index_inf;         // for a valve of index k, indexing this matrix gives the x-index of the cell where this valve will create a source term in the hard-coded source domain. This is for the bottom wall.
	int  *mfr_index_sup;        // for a valve of index k, indexing this matrix gives the x-index of the cell where this valve will create a source term in the hard-coded source domain. This is for the top wall.
	int  *mfr_n;                // the amount of source cells that this valve has, between the start- and end.
    
    int *fem_index_inf;			// The index in the x-direction where a valve starts
	int  *fem_index_sup;		// The index in the x-direction where a valve ends. 
	int  *fem_n;                // The amount of cells between the ending- and starting associated pressure index.
    int **p_neighbour;          // For each read valve, and for every FEM node, x-index of fluid cell associated
    double **p_coef;            // Interpolation coefficients for pressure at FEM nodes. Unknown what it physically represents. populated in build_fem_interface()
    double *mfr;
	double *mean_p_sup;			// The average pressure in the region above a valve (in the 'top' domain);
	double *mean_p_inf;			// The average pressure in the region where the valve moves around in (in the 'low' domain).
	double *mean_rho_sup;		// The average density in the region above a valve (in the 'top' domain);
	double *ytip;				// The y-coordinate of the tip for each valve. After initialisation, is set using u2.
	double *pratio;				// The ratio of the average pressure in the sampling region above and inside the lower domain.
	double *stage_mfr;			// Mass flow rate through each reed valve
    int n_cell_p_fem;			// The amount of cells in the y-direction that are used to determine the mean/average value.
    char valvename[40], index_char[10];
    double mfr_tot = 0.0;

    /* ROUTINE TO PERFORM IF SOLID MODELLING IS ACTIVATED */
    if (SOLID_ON==1)
    {
        /* NUMBER OF CELLS TO CONSIDER ALONG Y FOR ESTIMATION OF PRESSURE (PRESSURE RATIO AND FEM MODEL) */
        n_cell_p_fem = (int)fmax(1.0,floor(L_T/5.0/(L_TUBE/NX)));
        printf("Pressure ratios based on %i cells.\n", n_cell_p_fem);

    	/* MEMORY ALLOCATION */
    	N_NODE = N_FEM+1;
    	N_DOF= N_DOF_PER_NODE*N_NODE;
    	x_FEM = init_matrix(N_VALVE,N_FEM,1);               // Just creates a matrix of N_NODE * N_VALVE * 1 to store the x locations of the valves.
    	y_FEM = init_matrix(N_VALVE,N_FEM,1);
    	K = init_matrix(N_DOF,N_DOF,0);                     
    	C = init_matrix(N_DOF,N_DOF,0);                     
    	M = init_matrix(N_DOF,N_DOF,0);                     
        R1 = init_matrix(N_DOF,N_DOF,0);                    
        R2 = init_matrix(N_DOF,N_DOF,0);                    // Newmark matrix 2
        R3 = init_matrix(N_DOF,N_DOF,0);                    // Newmark matrix 3
	    b = init_array(N_FEM,1);
	    h = init_array(N_FEM,1);
    	A = init_array(N_FEM,1);
    	I = init_array(N_FEM,1);
    	U0_DOF = init_matrix(N_VALVE,N_DOF,0);
        U1_DOF = init_matrix(N_VALVE,N_DOF,0);
        U2_DOF = init_matrix(N_VALVE,N_DOF,0);
        U2_DOF_K = init_matrix(N_VALVE,N_DOF,0);
    	p_FEM = init_matrix(N_VALVE,N_NODE,0);
    	F_DOF = init_matrix(N_VALVE,N_DOF,0);
    	p_coef = init_matrix(N_VALVE,N_NODE,0);
    	p_neighbour = init_int_matrix(N_VALVE,N_NODE,0);
    	mfr_index_inf = malloc(N_VALVE*sizeof(int));
    	mfr_index_sup = malloc(N_VALVE*sizeof(int));
    	mfr_n = malloc(N_VALVE*sizeof(int));
    	fem_index_inf = malloc(N_VALVE*sizeof(int));
    	fem_index_sup = malloc(N_VALVE*sizeof(int));
    	fem_n = malloc(N_VALVE*sizeof(int));
        mfr = malloc(N_VALVE*sizeof(double));
        mean_p_inf = malloc(N_VALVE*sizeof(double));
        mean_p_sup = malloc(N_VALVE*sizeof(double));
        mean_rho_sup = malloc(N_VALVE*sizeof(double));
        ytip = malloc(N_VALVE*sizeof(double));
        pratio = malloc(N_VALVE*sizeof(double));
        stage_mfr = malloc(N_VALVE*sizeof(double));


        /* COMPUTE DEGREES OF FREEDOM AND FEM MATRICES */
	    act_DOF = compute_solid(N_VALVE,N_NODE,x_FEM,X_V_START,b,h,A,I,B0,B1,H0,H1,L_V,L_FIX,N_FIX,N_CLAMP,N_DOF_PER_NODE,&N_ACTIVE,&N_INACTIVE);
	    L_K = init_matrix(N_ACTIVE,N_ACTIVE,0);
	    LT_K = init_matrix(N_ACTIVE,N_ACTIVE,0);
	    L_R1 = init_matrix(N_ACTIVE,N_ACTIVE,0);
	    LT_R1 = init_matrix(N_ACTIVE,N_ACTIVE,0);
        printf("\nNumber of active/inactive FEM DOFs: %d / %d.\n", N_ACTIVE, N_INACTIVE);

        /* COMPUTE FEM MATRICES */
	    build_mass_mat(N_FEM,M,x_FEM[0],b,h,A,I,B0,B1,H0,H1,L_V,RHO_V,N_DOF_PER_NODE);
	    build_stiff_mat(N_FEM,K,x_FEM[0],b,h,A,I,B0,B1,H0,H1,L_V,E_V,N_DOF_PER_NODE);
	    build_damp_mat(N_DOF,C,K,M,RAYLEIGH_ALPHA,RAYLEIGH_BETA);
        build_newmark_mat(N_DOF,C,K,M,DT,R1,R2,R3);

        /* CHOLESKY DECOMPOSITION OF MASS AND STIFFNESS MATRICES */
        cholesky_decomposition(N_ACTIVE,R1,LT_R1,act_DOF);
        cholesky_decomposition(N_ACTIVE,K,LT_K,act_DOF);
        transpose(N_ACTIVE,LT_K,L_K);
        transpose(N_ACTIVE,LT_R1,L_R1);

        /* EXPORT ALL MATRICES IN .DAT FORMAT */
        export_fem_data(N_FEM,N_DOF,x_FEM[0],b,h,A,I,act_DOF,K,M,C,FEM_SECTIONS,FEM_K,FEM_M,FEM_C,EXP_EXTENSION);
        export_valve_data(N_VALVE,N_FEM+1,x_FEM,y_FEM,R0,0.0,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT,NGHOST);
        printf("\nFEM matrices and FEM model computed and exported...\n");

	    /* INITIAL FEM MESH AND INDEXES NECESSARY TO APPLY LOADS AND SOURCE TERMS */
	    for (int k = 0; k < N_VALVE; ++k)
	    {
	    	// Locate cells where the source term must be applied
	    	mfr_index_inf[k] = 0;
	    	mfr_index_sup[k] = 0;
            
            // Sets the x-index where the source term cells due to the presence of a valve start. 
	    	while ( (mfr_index_inf[k] < NXtot[dom_low]-NGHOST) && (x[dom_low][mfr_index_inf[k]][NGHOST] < X_V_START[k]+L_FIX+L_HOLE*(1.0-HOLE_FACTOR)) )
	    	{
	    		mfr_index_inf[k]++;
	    	}
            // Sets the x-index where the source term cells due to the presence of a valve end.
	    	while (mfr_index_sup[k]<NXtot[dom_low]-NGHOST && x[dom_low][mfr_index_sup[k]+1][NGHOST]<X_V_START[k]+L_FIX+L_HOLE)
	    	{
	    		mfr_index_sup[k]++;
	    	}
            // the amount of source cells that this valve has, between the start- and end.
	    	mfr_n[k] = 1 + mfr_index_sup[k] - mfr_index_inf[k];

	    	// Locate cells where the pressure field can be used to compute forces on the reed petal.
            // A: I think it just sets some cells at the tip.
	    	fem_index_inf[k] = 0;
	    	fem_index_sup[k] = 0;
            // But why is it offset by 2?
	    	while ( fem_index_inf[k]<NXtot[dom_low]-NGHOST && x[dom_low][fem_index_inf[k]+2][NGHOST]<X_V_START[k] )
	    	{
	    		fem_index_inf[k]++;
	    	}
	    	while ( fem_index_sup[k]<NXtot[dom_low]-NGHOST && x[dom_low][fem_index_sup[k]+1][NGHOST]<X_V_START[k]+L_T )
	    	{
	    		fem_index_sup[k]++;
	    	}
	    	fem_n[k] = 1 + fem_index_sup[k] - fem_index_inf[k];

	    	// Associate to FEM nodes a couple of fluid cells based on which node pressure is interpolated
	    	build_fem_interface(N_NODE,fem_n[k],fem_index_inf[k],x[dom_low],xc[dom_low],x_FEM[k], p_neighbour[k], p_coef[k], NGHOST);

            // Estimate pressure ratio and density around the valves
            mean_p_sup[k] = mean_at_valve(fem_index_inf[k],fem_n[k],NGHOST,n_cell_p_fem,p[dom_up]);
            mean_rho_sup[k] = mean_at_valve(fem_index_inf[k],fem_n[k],NGHOST,n_cell_p_fem,rho[dom_up]);
            mean_p_inf[k] = mean_at_valve(fem_index_inf[k],fem_n[k],NYtot[dom_low]-NGHOST-n_cell_p_fem,n_cell_p_fem,p[dom_low]);
	    }

        // DEFINE INTERFACE PRESSURE DIFFERENCE VECTOR
        // dp_interface = init_interface(N_VALVE,fem_n);

        /* EXPORT INITIAL MASS FLOW RATE AND TIP DISPLACEMENT FOR ALL VALVES (SET TO 0) */
        for (int k = 0; k < N_VALVE; ++k)
        {
            // Compute pressure at FEM nodes and deduce load on each node
            // fem_pressure(fem_n[k],fem_index_inf[k],x_FEM[k],NYtot[dom_low],n_cell_p_fem,p[dom_low],p[dom_up],NGHOST,dp_interface[k]);
            // fem_load(N_FEM,N_DOF,N_DOF_PER_NODE,N_CLAMP,NGHOST,X_V_START[k],p_neighbour[k],fem_n[k],dp_interface[k],x[dom_low],x_FEM[k],B0,B1,L_T,F_DOF[k]);
            // damping_load(N_FEM,N_DOF,N_DOF_PER_NODE,RHO_V,L_V,H0,H1,B0,B1,F_DOF[k]);

        	// Setting p_FEM, the pressure differences between the up- and the bottom.
        	
            fem_pressure(N_NODE,fem_n[k],fem_index_inf[k],x_FEM[k],p_neighbour[k],p_coef[k],NYtot[dom_low],p[dom_low],p[dom_up],NGHOST,p_FEM[k]);
        	// Sets F_DOF for the current valve.
            fem_load(N_FEM,N_DOF,N_CLAMP,N_DOF_PER_NODE,p_FEM[k],U2_DOF[k],F_DOF[k],x_FEM[k],b);

            // Solve initial FEM problem based on static assumption with initial pressure field
            cholesky_solve(N_DOF,N_ACTIVE,L_K,F_DOF[k],act_DOF,U2_DOF[k]);
            update_valve(N_NODE,N_DOF_PER_NODE,U0_DOF[k],U1_DOF[k],U2_DOF[k],y_FEM[k]);
            for (int k = 0; k < N_VALVE; ++k)
            {
            	
                for (int i = 0; i < N_DOF; ++i)
                {
                	// Artificially set the solution for the previous time step to be equal to this one; so it starts in stand still.
                    U1_DOF[k][i] = U2_DOF[k][i];
                    U0_DOF[k][i] = U2_DOF[k][i];
                }

                for (int i = 0; i < N_DOF; ++i)
                {
                    F_DOF[k][i]=0.0;
                }
            }

            ytip[k] = U2_DOF[k][N_DOF_PER_NODE*(N_NODE-1)]; // Set based on y deflection for each valve.
            stage_mfr[k] = 0.0;
            pratio[k] = mean_p_inf[k]/mean_p_sup[k];

            // strcpy(valvename,Y_TIP_FILENAME);
            // sprintf(index_char,"_%i",k);
            // strcat(valvename,index_char);
            // append_data_to_file(p_wall,0.0,P_MEAN_WALL_FILENAME,EXP_EXTENSION);
            // append_data_to_file(U2_DOF[k][N_DOF_PER_NODE*(N_NODE-1)],0.0,valvename,EXP_EXTENSION);
            // strcpy(valvename,MFR_TOT_FILENAME);
            // sprintf(index_char,"_%i",k);
            // strcat(valvename,index_char);
            // append_data_to_file(0.0,0.0,valvename,EXP_EXTENSION);
            // append_data_to_file(0.0,0.0,MFR_TOT_FILENAME,EXP_EXTENSION);
            // append_data_to_file(0.0,0.0,PFR_FILENAME,EXP_EXTENSION);
            // strcpy(valvename,P_RATIO_FILENAME);
            // sprintf(index_char,"_%i",k);
            // strcat(valvename,index_char);
            // append_data_to_file(0.0,0.0,valvename,EXP_EXTENSION);    
            // append_data_to_file(plenum_p,0.0,PLEN_P_FILENAME,EXP_EXTENSION);
        }

        // Export tip displacement, pressure ratio and mass flow rate
        append_multidata_to_file(N_VALVE,ytip,0.0,Y_TIP_FILENAME,EXP_EXTENSION);
        append_multidata_to_file(N_VALVE,stage_mfr,0.0,MFR_FILENAME,EXP_EXTENSION);
        append_multidata_to_file(N_VALVE,pratio,0.0,P_RATIO_FILENAME,EXP_EXTENSION);

        printf("\nInitial solid conditions applied...\n");
    }


    /* COUNTING TIME BETWEEN ITERATIONS */
    clock_gettime(CLOCK_REALTIME,&start_time_loop);


    /* TIME LOOP VARIABLES */
    double CFL=0;                        // For display of CFL
    double SIM_TIME_LOOP;                // To monitor real time between two exports
    double total_mfr = 0.0, total_pfr = 0.0;

    /* CELL VOLUME AT REED VALVE INTERFACE BETWEEN INSIDE/OUTSIDE TUBE */
    // double v_cell = M_PI*(pow(y[dom_low][NGHOST][NYtot[dom_low]-NGHOST-1],2)-pow(y[dom_low][NGHOST][NYtot[dom_low]-NGHOST-2],2))*(x[dom_low][NGHOST+1][NYtot[dom_low]-NGHOST-2]-x[dom_low][NGHOST][NYtot[dom_low]-NGHOST-1]);
    double r_low = y[dom_low][NGHOST][NYtot[dom_low]-NGHOST]; // I think this tracks the height of the bottom wall.
    double r_up = y[dom_up][NGHOST][NGHOST]; // I think this tracks the height of the upper bottom wall.
    // printf("Interface cell volume is: %g m^3.\n",v_cell);


    printf("Starting main time loop...\n");
    /* START TIME LOOP */
    for (int currentTimeStep = 1; currentTimeStep < Ntstep; ++currentTimeStep)
    {
        /* 4TH ORDER RUNGE-KUTTA PREDICTOR-CORRECTOR LOOP */
        for (int currentRungeKuttaIter = 0; currentRungeKuttaIter < RK_ORDER; ++currentRungeKuttaIter)
        {
            /* FOR FIRST ITERATION, TAKE PARAMETERS OF PREVIOUS SOLUTION */
            for (int k = 0; k < NDOMAIN; ++k)
            {
                #pragma omp parallel for
                for (int i = 0; i < NXtot[k]; ++i)
                {
                    for (int j = 0; j < NYtot[k]; ++j)
                    {
                    	
                        // Take solution from previous iteration for this new iteration
                        if(currentRungeKuttaIter==0)
                        {
                            // First Runge-Kutta iteration: Copy previous time step into Runge-kutta buffer.
                            rhoRK[k][i][j]=rho[k][i][j];
                            uRK[k][i][j]=u[k][i][j];
                            vRK[k][i][j]=v[k][i][j];
                            ERK[k][i][j]=E[k][i][j];
                            pRK[k][i][j]=p[k][i][j];
                            TRK[k][i][j]=T[k][i][j];
                            HRK[k][i][j]=H[k][i][j];
                        }
                        else
                        {
                            // Second, third, etc. iterations...
                        	// Copy 'new solution buffer' into runge kutta buffer.
                            rhoRK[k][i][j]=rhot[k][i][j];
                            uRK[k][i][j]=ut[k][i][j];
                            vRK[k][i][j]=vt[k][i][j];
                            ERK[k][i][j]=Et[k][i][j];
                            pRK[k][i][j]=pt[k][i][j];
                            TRK[k][i][j]=Tt[k][i][j];
                            HRK[k][i][j]=Ht[k][i][j];
                        }

                        // Reset all variable arrays used for final memory saving
                        rhot[k][i][j]=0;
                        ut[k][i][j]=0;
                        vt[k][i][j]=0;
                        Et[k][i][j]=0;
                        pt[k][i][j]=0;
                        Tt[k][i][j]=0;
                        Ht[k][i][j]=0;

                        // In addition, reset "sonicPoints" arrays
                        sonic_x[k][i][j]=0;
                        sonic_y[k][i][j]=0;
                    }
                }
            }

            /* SOLVE SOLID AND COMPUTE SOURCE TERM BASED ON CURRENT FLUID SOLUTION */
            if (SOLID_ON==1)
            {
                #pragma omp parallel for num_threads(N_VALVE)
                for (int k = 0; k < N_VALVE; ++k)
                {
                    // Compute pressure difference at FEM nodes and deduce load on each element/node
                    // fem_pressure(fem_n[k],fem_index_inf[k],x_FEM[k],NYtot[dom_low],n_cell_p_fem,pRK[dom_low],pRK[dom_up],NGHOST,dp_interface[k]);
                    // fem_load(N_FEM,N_DOF,N_DOF_PER_NODE,N_CLAMP,NGHOST,X_V_START[k],p_neighbour[k],fem_n[k],dp_interface[k],x[dom_low],x_FEM[k],B0,B1,L_T,F_DOF[k]);

                    fem_pressure(N_NODE,fem_n[k],fem_index_inf[k],x_FEM[k],p_neighbour[k],p_coef[k],NYtot[dom_low],pRK[dom_low],pRK[dom_up],NGHOST,p_FEM[k]);
                	// Populate F_DOF 
                    fem_load(N_FEM,N_DOF,N_CLAMP,N_DOF_PER_NODE,p_FEM[k],U2_DOF[k],F_DOF[k],x_FEM[k],b);
                    fem_flow_damping(N_FEM,N_DOF,N_CLAMP,N_DOF_PER_NODE,U1_DOF[k],U2_DOF[k],b,h,RHO_V,F0,DT,C1,C2,C3,F_DOF[k]);

                    // Solve FEM system and obtain Runge-Kutta valve distorsion at current RK loop
                    newmark_solve(N_DOF,N_ACTIVE,L_R1,R2,R3,F_DOF[k],act_DOF,U1_DOF[k],U2_DOF[k],U2_DOF_K[k]);

                    // Update FEM mesh and contrain displacement in positive domain
                    update_valve(N_NODE,N_DOF_PER_NODE,U0_DOF[k],U1_DOF[k],U2_DOF_K[k],y_FEM[k]);
                }
            }

            if (SOLID_ON==1)
            {
                // #pragma omp parallel for num_threads(N_VALVE)
                for (int k = 0; k < N_VALVE; ++k)
                {
                    // Compute mean pressure on each side and mass-flow rate at valve
                    mean_p_inf[k] = mean_at_valve(fem_index_inf[k],fem_n[k],NYtot[dom_low]-NGHOST-n_cell_p_fem,n_cell_p_fem,pRK[dom_low]);
                    mean_p_sup[k] = mean_at_valve(fem_index_inf[k],fem_n[k],NGHOST,n_cell_p_fem,pRK[dom_up]);
                    mean_rho_sup[k] = mean_at_valve(fem_index_inf[k],fem_n[k],NGHOST,n_cell_p_fem,rhoRK[dom_up]);
                    // double cur_p_inf = mean_at_valve(fem_index_inf[k],fem_n[k],NYtot[dom_low]-NGHOST-n_cell_p_fem,n_cell_p_fem,pRK[dom_low]);
                    // double cur_p_sup = mean_at_valve(fem_index_inf[k],fem_n[k],NGHOST,n_cell_p_fem,pRK[dom_up]);
                    // double cur_rho_sup = mean_at_valve(fem_index_inf[k],fem_n[k],NGHOST,n_cell_p_fem,rhoRK[dom_up]);

                    // Tip deflection and angle
                    double y_tip = y_FEM[k][N_FEM-1];

                    if (isnan(mean_p_sup[k]) || isnan(mean_p_inf[k]) || isnan(mean_rho_sup[k]))
                    {
                        printf("VALVE %d: PINF=%f, PSUP=%f, RATIO=%f\n",k,mean_p_inf[k],mean_p_sup[k],mean_p_inf[k]/mean_p_sup[k]);
                        printf("Crash at t=%f, ytip=%f\n",currentTimeStep*DT,y_tip);
                    }

                    // Time-average value of mean fields and mass flow rate
                	// A: SET MASS FLOW RATE, STORED for use in compute_cell_source
                    mfr[k] = compute_mfr(mean_p_inf[k],mean_p_sup[k],mean_rho_sup[k],y_tip,L_T,B0,B1,GAMMA,R);
                    
                }
            }

            /* GENERATE VALUES IN GHOST CELLS BASED ON PREVIOUS ITERATION */
            update_ghost_cells(NDOMAIN,NXtot,NYtot,x,y,xold,yold,rhoRK,uRK,vRK,pRK,HRK,ERK,TRK,B_LOC,B_TYPE,R0,NGHOST,GAMMA,R,P0,T0,M0); // As the runge-kutta buffer is apparently the one we're storing all the information in, this is the one we're populating.
            
            /* IDENTIFY SHOCK FRONTS IN EACH DIRECTION BEFORE PERFORMING AUSM-DV FLUX SPLITTING */
            for (int k = 0; k < NDOMAIN; ++k)
            {
                #pragma omp parallel for
                for (int i = NGHOST; i < NXtot[k]-NGHOST; ++i)
                {
                    double rhoL,uL,vL,pL,rhoR,uR,vR,pR,cL,cR;

                    for (int j = NGHOST; j < NYtot[k]-NGHOST; ++j)
                    {
                        rhoL=0.0;uL=0.0;vL=0.0;pL=0.0;rhoR=0.0;uR=0.0;vR=0.0;pR=0.0;cL=0.0;cR=0.0;
                    	
                        // MUSCL ON RIGHT FACE -> Right face flux
                        rhoL = MUSCL(rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],rhoRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],uRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],pRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],rhoRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],uRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],pRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);

                    	// If on this face, one normal direction is supersonic and the other is not, set the flag.
                        if( ( (uL-cL)>0 && (uR-cR)<0 ) || ( (uL+cL)>0 && (uR+cR)<0 ) )
                        {
                            // A sonic point is detected in the X direction between i and i+1
                            #pragma omp atomic
                            sonic_x[k][i][j] += 1;
                            #pragma omp atomic
                            sonic_x[k][i+1][j] += 1;
                        }

                        // MUSCL ON LEFT FACE -> Left face flux
                        rhoL = MUSCL(rhoRK[k][i-2][j],rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[k][i-2][j],uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i-2][j],pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i-2][j],rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[k][i-2][j],uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i-2][j],pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);
                        if( ( (uL-cL)>0 && uR-cR<0 ) || ( uL+cL>0 && uR+cR<0 ) )
                        {
                            // A sonic point is detected between i-1 and i
                            #pragma omp atomic
                            sonic_x[k][i-1][j] += 1;
                            #pragma omp atomic
                            sonic_x[k][i][j] += 1;
                        }

                        // MUSCL ON TOP FACE -> Top face flux (Left = Down and Right = Up)
                        rhoL = MUSCL(rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],rhoRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],vRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],pRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],rhoRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],vRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],pRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);
                        if( ( (vL-cL)>0 && (vR-cR)<0 ) || ( (vL+cL)>0 && (vR+cR)<0 ) )
                        {
                            // A sonic point is detected between j and j+1
                            #pragma omp atomic
                            sonic_y[k][i][j] += 1;
                            #pragma omp atomic
                            sonic_y[k][i][j+1] += 1;
                        }
                        // MUSCL ON DOWN FACE -> Down face flux (Left = Down and Right = Up)
                        rhoL = MUSCL(rhoRK[k][i][j-2],rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[k][i][j-2],vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i][j-2],pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i][j-2],rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[k][i][j-2],vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i][j-2],pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);
                        if( ( (vL-cL)>0 && (vR-cR)<0 ) || ( (vL+cL)>0 && (vR+cR)<0 ) )
                        {
                            // A sonic point is detected between j-1 and j
                            #pragma omp atomic
                            sonic_y[k][i][j-1] += 1;
                            #pragma omp atomic
                            sonic_y[k][i][j] += 1;
                        }
                    }
                }
            }


            /* SOLVING ALL DOMAINS WITH AUSM-DV OR HANEL'S SCHEME DEPENDING ON THE POSITION OF SHOCK FRONTS */
            for (int k = 0; k < NDOMAIN; ++k)
            {
                #pragma omp parallel for
                for (int i = NGHOST; i < NXtot[k]-NGHOST; ++i)
                {
                    double RK_k;
                    
                	
					double uL; // MUSCL of u-velocity, left face.
					double vL; // MUSCL of v-velocity, left face.
					double pL; // MUSCL of pressure, left face.
                	double rhoL; // MUSCL of density, left face.
                	double EL; // MUSCL of energy, left face.
                	double HL; // MUSCL of enthalpy(?), left face.
					
					double uR; // MUSCL of u-velocity, right face.
					double vR; // MUSCL of v-velocity, right face.
					double pR; // MUSCL of pressure, right face.
                	double rhoR; // MUSCL of density, right face.
                	double ER; // MUSCL of energy, right face.
                	double HR; // MUSCL of enthalpy(?), right face.
					
					double  theta;
                	
                    double flux_l[4]; // flux terms of over the left face in the current cell index. Index corresponds to euler equation: [0 : density term, 1: x-axis velocity, 2: y/r-axis velocity, 3: specific internal energy]
					double flux_r[4]; // flux terms of over the right face in the current cell index. Index corresponds to euler equation: [0 : density term, 1: x-axis velocity, 2: y/r-axis velocity, 3: specific internal energy]
					double flux_u[4]; // flux terms of over the up face in the current cell index. Index corresponds to euler equation: [0 : density term, 1: x-axis velocity, 2: y/r-axis velocity, 3: specific internal energy]
					double flux_d[4]; // flux terms of over the down face in the current cell index. Index corresponds to euler equation: [0 : density term, 1: x-axis velocity, 2: y/r-axis velocity, 3: specific internal energy]
					double flux[4];
                	
                    double dx; // The distance between this cell center and the one at index x+1
					double dr; // The distance between this cell center and the one at index x+1
                	
                    double y_tip; // deflection of the tip of the reed valve
					double  theta_tip; // angle that the tip makes to its x-axis.
					double  source[4];
					double  dx_source;
					double  dy_source;
                	
                    int m;

                    for (int j = NGHOST; j < NYtot[k]-NGHOST; ++j)
                    {   
                        RK_k=0.0;
                        rhoL=0.0;uL=0.0;vL=0.0;pL=0.0;rhoR=0.0;uR=0.0;vR=0.0;pR=0.0;EL=0.0;ER=0.0;HL=0.0;HR=0.0, theta=0.0;
                        for (int p = 0; p < 4; ++p)
                        {
                            flux_l[p]=0;
                            flux_r[p]=0;
                            flux_u[p]=0;
                            flux_d[p]=0;
                            source[p]=0;
                            flux[p]=0;
                        }
                        dx=0.0;dr=0.0;
                        y_tip=0.0;
                        theta_tip=0.0;
                        m=0;

                    	// Calculation procedure MUSCL fluxes is identical to the one above, but now it gets vL as well as uL. so in 2D.

                        // MUSCL ON RIGHT FACE -> Right face flux
                        rhoL = MUSCL(rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],rhoRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],uRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[k][i-1][j],vRK[k][i][j],vRK[k][i+1][j],vRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],pRK[k][i+2][j],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],rhoRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],uRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[k][i-1][j],vRK[k][i][j],vRK[k][i+1][j],vRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],pRK[k][i+2][j],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                    	// A: This looks unphysical, as the value can be set by a flux in a different direction, but still cause hanel to be used.
                        // FLUX SPLITTING ON RIGHT FACE based on sonic points in Y-direction
                        if (sonic_y[k][i][j]+sonic_y[k][i+1][j]>=0) // If either of the two sampled points has the 'sonic' flag set.
                        {
                            HANEL(flux_r,'H',rhoL,rhoR,uL,uR,vL,vR,pL,pR,HL,HR,R,GAMMA,ENTRO_FIX_C);
                        }
                        else
                        {
                            // AUSM-DV SCHEME FOR RIGHT FACE
                            AUSM_DV(flux_r,'H',rhoL,rhoR,uL,uR,vL,vR,pL,pR,HL,HR,R,GAMMA,AUSM_K,ENTRO_FIX_C);
                        }

                        rhoL=0,uL=0,vL=0,pL=0,HL=0;
                        rhoR=0,uR=0,vR=0,pR=0,HR=0;

                    	// These are all offset 1 cell to the left.
                        // MUSCL ON LEFT FACE -> Left face flux
                        rhoL = MUSCL(rhoRK[k][i-2][j],rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[k][i-2][j],uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[k][i-2][j],vRK[k][i-1][j],vRK[k][i][j],vRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i-2][j],pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i-2][j],rhoRK[k][i-1][j],rhoRK[k][i][j],rhoRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[k][i-2][j],uRK[k][i-1][j],uRK[k][i][j],uRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[k][i-2][j],vRK[k][i-1][j],vRK[k][i][j],vRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i-2][j],pRK[k][i-1][j],pRK[k][i][j],pRK[k][i+1][j],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // FLUX SPLITTING ON LEFT FACE based on sonic points in Y-direction
                        if (sonic_y[k][i][j]+sonic_y[k][i-1][j]>=0)
                        {
                            HANEL(flux_l,'H',rhoL,rhoR,uL,uR,vL,vR,pL,pR,HL,HR,R,GAMMA,ENTRO_FIX_C);
                        }
                        else
                        {
                            // AUSM-DV SCHEME FOR LEFT FACE
                            AUSM_DV(flux_l,'H',rhoL,rhoR,uL,uR,vL,vR,pL,pR,HL,HR,R,GAMMA,AUSM_K,ENTRO_FIX_C);
                        }

                        rhoL=0,uL=0,vL=0,pL=0,HL=0;
                        rhoR=0,uR=0,vR=0,pR=0,HR=0;
                        theta = 0.0;

                        // VERTICAL FLUX-VECTOR SPLITTING: G(1/2) and G(-1/2)
                    	// A: Literally same, but now the last index, aka y.
                        
                        // MUSCL ON TOP FACE -> Top face flux (Left = Down and Right = Up)
                        rhoL = MUSCL(rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],rhoRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[k][i][j-1],uRK[k][i][j],uRK[k][i][j+1],uRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],vRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],pRK[k][i][j+2],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],rhoRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[k][i][j-1],uRK[k][i][j],uRK[k][i][j+1],uRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],vRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],pRK[k][i][j+2],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // PROJECT VELOCITIES IN THE FRAME OF THE UPPER FACE
                        // if (strcmp(SIM_CASE,"sup_plen_rocket")==0 || strcmp(SIM_CASE,"sup_plenum")==0)
                        // {
                        //     if (k==cone_index && i==NGHOST)
                        //     {
                        //         theta = atan((y[k][i][j]-y[k][i-1][j])/(x[k][i][j]-x[k][i-1][j]));
                        //         theta_frame_velocity(uL,vL,theta);
                        //         theta_frame_velocity(uR,vR,theta);
                        //     }
                        // }

                        // FLUX SPLITTING ON TOP FACE based on sonic points in X-direction
                        if (sonic_x[k][i][j]+sonic_x[k][i][j+1]>=0)
                        {
                            // HANEL SCHEME FOR TOP FACE
                        	// Note that u & v are flipped, in addition to the 'v' flag being applied.
                            HANEL(flux_u,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,ENTRO_FIX_C);
                        }
                        else
                        {
                        	// Note that u & v are flipped, in addition to the 'v' flag being applied.
                            // AUSM-DV SCHEME FOR TOP FACE
                            AUSM_DV(flux_u,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,AUSM_K,ENTRO_FIX_C);
                        }

                        rhoL=0,uL=0,vL=0,pL=0,HL=0;
                        rhoR=0,uR=0,vR=0,pR=0,HR=0;
                        theta = 0.0;

                    	// Again, literally offset by -1
                        // MUSCL ON DOWN FACE -> Down face flux (Left = Down and Right = Up)
                        rhoL = MUSCL(rhoRK[k][i][j-2],rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[k][i][j-2],uRK[k][i][j-1],uRK[k][i][j],uRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[k][i][j-2],vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[k][i][j-2],pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[k][i][j-2],rhoRK[k][i][j-1],rhoRK[k][i][j],rhoRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[k][i][j-2],uRK[k][i][j-1],uRK[k][i][j],uRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[k][i][j-2],vRK[k][i][j-1],vRK[k][i][j],vRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[k][i][j-2],pRK[k][i][j-1],pRK[k][i][j],pRK[k][i][j+1],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // PROJECT VELOCITIES IN THE FRAME OF THE CURRENT FACE
                        // if (strcmp(SIM_CASE,"sup_plen_rocket")==0 || strcmp(SIM_CASE,"sup_plenum")==0)
                        // {
                        //     if (k==cone_index && i==NGHOST)
                        //     {
                        //         theta = atan((y[k][i][j-1]-y[k][i-1][j-1])/(x[k][i][j-1]-x[k][i-1][j-1]));
                        //         theta_frame_velocity(uL,vL,theta);
                        //         theta_frame_velocity(uR,vR,theta);
                        //     }
                        // }

                        // FLUX SPLITTING ON DOWN FACE based on sonic points in X-direction
                        if (sonic_y[k][i][j]+sonic_y[k][i][j-1]>=0)
                        {
                        	// Note that u & v are flipped, in addition to the 'v' flag being applied.
                            // HANEL SCHEME FOR DOWN FACE
                            HANEL(flux_d,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,ENTRO_FIX_C);
                        }
                        else
                        {
                            // AUSM-DV SCHEME FOR DOWN FACE
                        	// Note that u & v are flipped, in addition to the 'v' flag being applied.
                            AUSM_DV(flux_d,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,AUSM_K,ENTRO_FIX_C);
                        }

                        rhoL=0,uL=0,vL=0,pL=0,HL=0;
                        rhoR=0,uR=0,vR=0,pR=0,HR=0;

                        
                        /* COMPUTATION OF TOTAL FLUX - FINITE VOLUME FORMULATION */
                        // if (strcmp(SIM_CASE,"sup_plen_rocket")==0 || strcmp(SIM_CASE,"sup_plenum")==0)
                        // {
                        //     if (k==cone_index && i==NGHOST)
                        //     {
                        //         theta = atan((y[k][i][j]-y[k][i-1][j])/(x[k][i][j]-x[k][i-1][j]));
                        //         theta_frame_velocity(uL,vL,theta);
                        //         theta_frame_velocity(uR,vR,theta);
                        //     }
                        // }
                        // else
                        // {
                            dx = x[k][i+1][j] - x[k][i][j];
                            dr = y[k][i][j+1] - y[k][i][j];

                    		// Accumulation because of what goes in/out
                            flux[0] = ((flux_r[0] - flux_l[0])/dx + (flux_u[0] - flux_d[0])/dr)*DT;
                            flux[1] = ((flux_r[1] - flux_l[1])/dx + (flux_u[1] - flux_d[1])/dr)*DT;
                            flux[2] = ((flux_r[2] - flux_l[2])/dx + (flux_u[2] - flux_d[2])/dr)*DT;
                            flux[3] = ((flux_r[3] - flux_l[3])/dx + (flux_u[3] - flux_d[3])/dr)*DT;
                        // }

                        /* SOURCE TERM CORRESPONDING TO CYLINDRICAL SYSTEM OF COORDINATES */
                    	// A: This is the Hr
                        flux[0] += rho[k][i][j]*v[k][i][j]/yc[k][i][j]*DT;
                        flux[1] += rho[k][i][j]*v[k][i][j]*u[k][i][j]/yc[k][i][j]*DT;
                        flux[2] += rho[k][i][j]*v[k][i][j]*v[k][i][j]/yc[k][i][j]*DT;
                        flux[3] += rho[k][i][j]*v[k][i][j]*(H[k][i][j])/yc[k][i][j]*DT;

                        /* IF NECESSARY, ADD/REMOVE FLUX AT INTERFACES */
                        if (SOLID_ON==1 && ((k==dom_low && j==NYtot[dom_low]-NGHOST-1) || (k==dom_up && j==NGHOST)))
                        {
                            // Locate associated valve and source term value
                            m = 0;
                            while( !(i>= mfr_index_inf[m] && i<=mfr_index_sup[m]) && m<N_VALVE)
                            {
                                m++;
                            }

                            // Add or remove according to domain (inverse convention because eventually substracted)
                            if (m<N_VALVE)
                            {
                                dx_source = x[k][i+1][j]-x[k][i][j];
                                dy_source = y[k][i][j+1]-y[k][i][j];
                                y_tip = U2_DOF[m][N_DOF_PER_NODE*(N_NODE-1)];
                                theta_tip = U2_DOF_K[m][N_DOF_PER_NODE*(N_NODE-1)+1];
                                compute_cell_source(k,dom_low,dom_up,r_low,r_up,mfr[m],y_tip,theta_tip,dy_source, N_V_PER_STAGE,GAMMA,R0,L_HOLE,HOLE_FACTOR,dx_source,B_HOLE,mfr_n[m],mean_rho_sup[m],mean_p_sup[m],mean_p_inf[m],rhoRK[dom_up][i][j],uRK[dom_up][i][j],vRK[dom_up][i][j],pRK[dom_up][i][j],B0,B1,L_T,source);

                            	// Extra source due to valve opening.
                                flux[0] -= source[0]*DT;
                                flux[1] -= source[1]*DT;
                                flux[2] -= source[2]*DT;
                                flux[3] -= source[3]*DT;
                                // printf("%f %f %f %f \n",source[0],source[1],source[2],source[3]);

                                if (isnan(source[0]) || isnan(source[1]) || isnan(source[2]) || isnan(source[3]))
                                {
                                    printf("DOM UP IS %d\n", dom_up);
                                    fprintf(stderr,"\nERROR:SOURCE IS NaN FOR VALVE %d at time %d and at domain %d %d %d: %f %f %f %f (ytip is %f)!\n",m,currentTimeStep,k,i,j,source[0],source[1],source[2],source[3],y_tip);
                                    exit(0);
                                }
                            }
                        }

                        /* INCREMENT RUNGE-KUTTA VARIABLES WITH CORRESPONDING FLUX */
                        RK_k = 1./(RK_ORDER-currentRungeKuttaIter);
                        rhot[k][i][j] = rho[k][i][j] - RK_k*flux[0];
                        ut[k][i][j] = (rho[k][i][j]*u[k][i][j] - RK_k*flux[1])/rhot[k][i][j];
                        vt[k][i][j] = (rho[k][i][j]*v[k][i][j] - RK_k*flux[2])/rhot[k][i][j];
                        Et[k][i][j] = E[k][i][j] - RK_k*flux[3];

                        /* DEDUCE REDUNDANT VARIABLES */
                        pt[k][i][j] = (GAMMA-1)*(Et[k][i][j] - 0.5*rhot[k][i][j]*(pow(ut[k][i][j],2)+pow(vt[k][i][j],2)));
                        Tt[k][i][j] = pt[k][i][j]/(R*rhot[k][i][j]);
                        Ht[k][i][j] = (Et[k][i][j]+pt[k][i][j])/rhot[k][i][j];
                    }
                }
            }
        } 
        /* END OF THE RUNGE-KUTTA LOOP */

        /* UPDATE SOLUTION BEFORE NEXT ITERATION */
        for (int k = 0; k < NDOMAIN; ++k)
        {
            #pragma omp parallel for num_threads(4) schedule(static)
            for (int i = 0; i < NXtot[k]; ++i)
            {
                #pragma omp parallel for num_threads(3) schedule(static)
                for (int j = 0; j < NYtot[k]; ++j)
                {
                    rho[k][i][j]=rhot[k][i][j];
                    u[k][i][j]=ut[k][i][j];
                    v[k][i][j]=vt[k][i][j];
                    E[k][i][j]=Et[k][i][j];
                    p[k][i][j]=pt[k][i][j];
                    T[k][i][j]=Tt[k][i][j];
                    H[k][i][j]=Ht[k][i][j];

                    if (isnan(p[k][i][j]))
                    {
                        fprintf(stderr,"SOLVED PRESSURE IS NaN at (%d,%d,%d)\n",k,i,j);
                    }
                }
            }
        }


        /* UPDATE SOLID NODES WITH LAST RUNGE-KUTTA SOLUTION */
        if (SOLID_ON==1)
        {
            for (int k = 0; k < N_VALVE; ++k)
            {
                for (int i = 0; i < N_DOF; ++i)
                {
                	// I think this shows that all U's are deflections at different time steps. U0 or U1
                    U0_DOF[k][i] = U1_DOF[k][i];
                    U1_DOF[k][i] = U2_DOF[k][i];
                    U2_DOF[k][i] = U2_DOF_K[k][i];
                }
            }
        }

        /* EXPORT SEVERAL FLUID-RELATED PARAMETERS */
        p_wall = 0.0, mfr_intake = 0.0, p_tube = 0.0, rho_tube = 0.0;
        // for (int i = NGHOST; i < NYtot[dom_low]-NGHOST; ++i)
        // {
        //     p_wall += p[dom_low][NGHOST][i]*(pow(y[dom_low][NGHOST][i+1],2)-pow(y[dom_low][NGHOST][i],2))/pow(R0,2);
        //     mfr_intake += rho[dom_low][NXtot[dom_low]-NGHOST-1][i]*u[dom_low][NXtot[dom_low]-NGHOST-1][i]*M_PI*(pow(y[dom_low][NXtot[dom_low]-NGHOST][i+1],2)-pow(y[dom_low][NXtot[dom_low]-NGHOST][i],2));
        // }
        p_wall = average_face(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],p[dom_low],NGHOST,NGHOST);
        mfr_intake = mfr_face(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],rho[dom_low],u[dom_low],NXtot[dom_low]-NGHOST-1,NGHOST);
        p_tube = domain_average(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],p[dom_low],NGHOST);
        rho_tube = domain_average(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],rho[dom_low],NGHOST);

        append_data_to_file(p_wall,DT*currentTimeStep,P_MEAN_WALL_FILENAME,EXP_EXTENSION);
        append_data_to_file(mfr_intake,DT*currentTimeStep,INTAKE_MFR_FILENAME,EXP_EXTENSION);
        append_data_to_file(p_tube,DT*currentTimeStep,P_TUBE_FILENAME,EXP_EXTENSION);
        append_data_to_file(rho_tube,DT*currentTimeStep,RHO_TUBE_FILENAME,EXP_EXTENSION);

        if (PLENUM_ON==1)
        {
            mfr_plenum = 0.0, p_plenum = 0.0, rho_plenum = 0.0, p_drag = 0.0;

            mfr_plenum = mfr_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],u[dom_up],NGHOST,NGHOST);
            p_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NGHOST);
            rho_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],NGHOST);
            p_drag = average_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NXtot[dom_up]-NGHOST-1,NGHOST);

            append_data_to_file(p_plenum,DT*currentTimeStep,PLEN_P_FILENAME,EXP_EXTENSION);
            append_data_to_file(rho_plenum,DT*currentTimeStep,PLEN_RHO_FILENAME,EXP_EXTENSION);
            append_data_to_file(mfr_plenum,DT*currentTimeStep,PLENUM_MFR_FILENAME,EXP_EXTENSION);
            append_data_to_file(p_drag,DT*currentTimeStep,PLENUM_DRAG_PRESSURE_FILENAME,EXP_EXTENSION);
        }


        /* EXPORT MASS FLOW RATE AND VALVE DISPLACEMENT IF REQUIRED */
        if (SOLID_ON==1)
        {
            // Total and local valve mass flow rate export
            total_mfr = 0.0;

            for (int i = 0; i < N_VALVE; ++i)
            {
                total_mfr += N_V_PER_STAGE*mfr[i];
                stage_mfr[i] = N_V_PER_STAGE*mfr[i];
                ytip[i] = U2_DOF[i][N_DOF_PER_NODE*(N_NODE-1)];
                pratio[i] = mean_p_inf[i]/mean_p_sup[i];
            }

            // PFR at current iteration
            if (strcmp(SIM_CASE,"sup_plen_rocket")!=0)
            {
                total_pfr += total_mfr*DT/(P0/R/T0*M_PI*R0*R0*L_TUBE);
            }
            else
            {
                total_pfr += total_mfr*DT/(p_prerun/R/T_prerun*M_PI*R0*R0*L_TUBE);
            }

            // PFR, stage mass flow rate and total mass flow rate export
            append_data_to_file(total_mfr,DT*currentTimeStep,MFR_TOT_FILENAME,EXP_EXTENSION);
            append_data_to_file(total_pfr,DT*currentTimeStep,PFR_FILENAME,EXP_EXTENSION);
            append_multidata_to_file(N_VALVE,stage_mfr,DT*currentTimeStep,MFR_FILENAME,EXP_EXTENSION);

            // Valve displacement and pressure ratio export
            append_multidata_to_file(N_VALVE,ytip,DT*currentTimeStep,Y_TIP_FILENAME,EXP_EXTENSION);
            append_multidata_to_file(N_VALVE,pratio,DT*currentTimeStep,P_RATIO_FILENAME,EXP_EXTENSION);
        }

        /* EXPORT ALL FLUID AND SOLID DATA EVERY N STEPS */
        if(currentTimeStep%N_EXPORT==0)
        {
            // Fluid data export
            export_fluid_data(NDOMAIN,NXtot,NYtot,x,y,xc,yc,rho,u,v,p,E,T,H,DT*currentTimeStep,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT);

            // Valve data export
            if (SOLID_ON==1)
            {
                // Valve displacement
            	export_valve_data(N_VALVE,N_FEM+1,x_FEM,y_FEM,R0,DT*currentTimeStep,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT,NGHOST);
            }
        }


        /* DISPLAY CFL */
        if(currentTimeStep%N_CFL==0)
        {
            // COMPUTE CFL (BASED ON TOTAL VELOCITY)
            CFL = get_cfl(x,y,u,v,T,R,GAMMA,NDOMAIN,NXtot,NYtot,DT);

            // COUNTING TIME BETWEEN ITERATIONS
            clock_gettime(CLOCK_REALTIME,&end_time_loop);
            SIM_TIME_LOOP = (end_time_loop.tv_sec - start_time_loop.tv_sec)+(end_time_loop.tv_nsec - start_time_loop.tv_nsec)/1E9;
            printf("Iteration %d (%.1f %%)... Max. CFL at t=%.5f sec is: %f. \n[completed in %f sec]\n",currentTimeStep,(double)(currentTimeStep)/Ntstep*100.0,DT*currentTimeStep,CFL,SIM_TIME_LOOP);
            clock_gettime(CLOCK_REALTIME,&start_time_loop);

            // RESET BEFORE NEXT LOOP
            SIM_TIME_LOOP = 0;
        }


    } // END OF TIME LOOP


    /* FINAL INFORMATION DISPLAY AT END OF SIMULATION */
    clock_gettime(CLOCK_REALTIME,&end_time);
    double SIM_TIME = (end_time.tv_sec - start_time.tv_sec)+(end_time.tv_nsec - start_time.tv_nsec)/1E9;
    printf("\nSimulation time was: %f s.",SIM_TIME);
    printf("\nSimulation time was: %f min.",SIM_TIME/60.0);
    printf("\nSimulation time was: %f hours.",SIM_TIME/3600.0);
    printf("\nSimulation time per cell and 1000 iterations was: %f s.",SIM_TIME/Ncell/Ntstep*1000);
    printf("\nSIMULATION COMPLETE! -Florian   ");

}


/* END OF FILE */

