/* 
* Rewrite of Florian's code, in cpp
*/

#include <ctime>
#include <iostream>

#include "colors.h"      // Colors for use in terminal
#include "parameters.h"  // List of parameters specified by user
#include "functions.h"   // General functions
#include "initial.h"     // Initial conditions suitable for each case
#include "muscl.h"       // MUSCL interpolation
#include "ausm.h"        // AUSM scheme
#include "microwave.h"   // Microwave heating and MSD theory
#include "export.h"      // Export functions
#include "bound_cond.h"  // Boundary conditions
#include "fem_model.h"   // FEM model of the reed valve

#include "sim_case.h"
#include "case_det_tube.cpp" // Temporary hard-code of specific case implementation



/* MAIN ROUTINE */
int main()
{
    /* STARTING PROGRAM, GETTING TIME... */
    std::cout << "Program starting...\n";
    std::time_t timeAtStartOfProgram = std::time(0);

    // INITIALIZE CASE
    SimCase simCase;
    LoadCase(&simCase);

    /* CREATE OUTPUT FOLDER */
    make_dir(OUT_FOLDERNAME);

    /* DISPLAY INFORMATION ON TERMINAL */
    std::cout << "Total length of the simulation is " << TSIM << "s, with dt = " << DT << " ms (" << simCase.totalSimulationTimeStepCount << " steps)." << std::endl;
    std::cout << "CFL will be displayed every " << simCase.runtimeParameters.numberOfIterationsBetweenCFLLog << "timesteps (" << simCase.runtimeParameters.numberOfIterationsBetweenCFLLog * DT*1000 <<"ms)." << std::endl;
    std::cout << "Field properties will be exported every " << simCase.runtimeParameters.numberOfIterationsBetweenDataExport << "timesteps (" << simCase.runtimeParameters.numberOfIterationsBetweenDataExport * DT*1000 <<"ms)." << std::endl;

    /* EXPORT PARAMETERS REQUIRED FOR POST-PROCESSING AND SOLUTION RECONSTRUCTION */
    //export_parameters(NDOMAIN,TSIM,DT,XSTART,YSTART,XLENGTH,YLENGTH,N_VALVE,N_FEM,NXtot,NYtot,NGHOST,N_EXPORT,W_FORMAT,PAR_FILENAME);


    /* INITIAL CONDITIONS BASED ON MICROWAVE DETONATION THEORY */
    ChapmanJougetDetonationSolution initialDetonationSolution = SolveChapmanJougetDetonationProblem(T0, P0, ETA, S0, R, GAMMA, L_TUBE, R0);
    std::cout << "Chapman-Jouget Detonation solution for initial conditions found after " << initialDetonationSolution.iters_performed << " iterations." << std::endl;


    /* INITIAL CONDITIONS ON DOMAINS */
    simCase.ApplyInitialConditions();
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
	int N_DOF, N_NODE, N_ACTIVE, N_INACTIVE, *act_DOF;
	double **x_FEM, **y_FEM, *b, *h, *A, *I;
    double **K, **C, **M, **L_K, **LT_K, **L_R1, **LT_R1, **R1, **R2, **R3;
    double **U0_DOF, **U1_DOF, **U2_DOF, **U2_DOF_K, **dp_interface, **p_FEM,**F_DOF;

    /* MASS-FLOW RATE RELATED PARAMETERS DUE TO VALVE */
    int *mfr_index_inf, *mfr_index_sup, *mfr_n; // To locate cells with source term
    int *fem_index_inf, *fem_index_sup, *fem_n; // To locate cells whose pressure is used
    int **p_neighbour; // For each FEM node, index of fluid cell associated (i and i+1)
    double **p_coef; // Interpolation coefficients for pressure at FEM nodes
    double *mfr, *mean_p_sup, *mean_p_inf, *mean_rho_sup, *ytip, *pratio, *stage_mfr;   // Mass flow rate through each reed valve
    int n_cell_p_fem;
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
    	x_FEM = init_matrix(N_VALVE,N_FEM,1);
    	y_FEM = init_matrix(N_VALVE,N_FEM,1);
    	K = init_matrix(N_DOF,N_DOF,0);
    	C = init_matrix(N_DOF,N_DOF,0);
    	M = init_matrix(N_DOF,N_DOF,0);
        R1 = init_matrix(N_DOF,N_DOF,0);
        R2 = init_matrix(N_DOF,N_DOF,0);
        R3 = init_matrix(N_DOF,N_DOF,0);
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
	    	while (mfr_index_inf[k]<NXtot[dom_low]-NGHOST && x[dom_low][mfr_index_inf[k]][NGHOST]<X_V_START[k]+L_FIX+L_HOLE*(1.0-HOLE_FACTOR))
	    	{
	    		mfr_index_inf[k]++;
	    	}
	    	while (mfr_index_sup[k]<NXtot[dom_low]-NGHOST && x[dom_low][mfr_index_sup[k]+1][NGHOST]<X_V_START[k]+L_FIX+L_HOLE)
	    	{
	    		mfr_index_sup[k]++;
	    	}
	    	mfr_n[k] = 1 + mfr_index_sup[k] - mfr_index_inf[k];

	    	// Locate cells where the pressure field can be used to compute forces on the reed petal
	    	fem_index_inf[k] = 0;
	    	fem_index_sup[k] = 0;
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

            fem_pressure(N_NODE,fem_n[k],fem_index_inf[k],x_FEM[k],p_neighbour[k],p_coef[k],NYtot[dom_low],p[dom_low],p[dom_up],NGHOST,p_FEM[k]);
            fem_load(N_FEM,N_DOF,N_CLAMP,N_DOF_PER_NODE,p_FEM[k],U2_DOF[k],F_DOF[k],x_FEM[k],b);

            // Solve initial FEM problem based on static assumption with initial pressure field
            cholesky_solve(N_DOF,N_ACTIVE,L_K,F_DOF[k],act_DOF,U2_DOF[k]);
            update_valve(N_NODE,N_DOF_PER_NODE,U0_DOF[k],U1_DOF[k],U2_DOF[k],y_FEM[k]);
            for (int k = 0; k < N_VALVE; ++k)
            {
                for (int i = 0; i < N_DOF; ++i)
                {
                    U1_DOF[k][i] = U2_DOF[k][i];
                    U0_DOF[k][i] = U2_DOF[k][i];
                }

                for (int i = 0; i < N_DOF; ++i)
                {
                    F_DOF[k][i]=0.0;
                }
            }

            ytip[k] = U2_DOF[k][N_DOF_PER_NODE*(N_NODE-1)];
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
    double r_low = y[dom_low][NGHOST][NYtot[dom_low]-NGHOST];
    double r_up = y[dom_up][NGHOST][NGHOST];
    // printf("Interface cell volume is: %g m^3.\n",v_cell);


    printf("Starting main time loop...\n");
    /* START TIME LOOP */
    for (int t = 1; t < Ntstep; ++t)
    {
        /* 4TH ORDER RUNGE-KUTTA PREDICTOR-CORRECTOR LOOP */
        for (int trk = 0; trk < RK_ORDER; ++trk)
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
                        if(trk==0)
                        {
                            // First Runge-Kutta iteration
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
                // #pragma omp parallel for num_threads(N_VALVE) (bug if activated?)
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
                        printf("Crash at t=%f, ytip=%f\n",t*DT,y_tip);
                    }

                    // Time-average value of mean fields and mass flow rate
                    mfr[k] = compute_mfr(mean_p_inf[k],mean_p_sup[k],mean_rho_sup[k],y_tip,L_T,B0,B1,GAMMA,R);
                    
                }
            }

            /* GENERATE VALUES IN GHOST CELLS BASED ON PREVIOUS ITERATION */
            update_ghost_cells(NDOMAIN,NXtot,NYtot,x,y,xold,yold,rhoRK,uRK,vRK,pRK,HRK,ERK,TRK,B_LOC,B_TYPE,R0,NGHOST,GAMMA,R,P0,T0,M0);
            
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
                    double rhoL,uL,vL,pL,rhoR,uR,vR,pR,EL,ER,HL,HR, theta;
                    double flux_l[4],flux_r[4],flux_u[4],flux_d[4],flux[4];
                    double dx,dr;
                    double y_tip, theta_tip, source[4], dx_source, dy_source;
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

                        // FLUX SPLITTING ON RIGHT FACE based on sonic points in Y-direction
                        if (sonic_y[k][i][j]+sonic_y[k][i+1][j]>=0)
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
                            HANEL(flux_u,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,ENTRO_FIX_C);
                        }
                        else
                        {
                            // AUSM-DV SCHEME FOR TOP FACE
                            AUSM_DV(flux_u,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,AUSM_K,ENTRO_FIX_C);
                        }

                        rhoL=0,uL=0,vL=0,pL=0,HL=0;
                        rhoR=0,uR=0,vR=0,pR=0,HR=0;
                        theta = 0.0;

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
                            // HANEL SCHEME FOR DOWN FACE
                            HANEL(flux_d,'V',rhoL,rhoR,vL,vR,uL,uR,pL,pR,HL,HR,R,GAMMA,ENTRO_FIX_C);
                        }
                        else
                        {
                            // AUSM-DV SCHEME FOR DOWN FACE
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

                            flux[0] = ((flux_r[0] - flux_l[0])/dx + (flux_u[0] - flux_d[0])/dr)*DT;
                            flux[1] = ((flux_r[1] - flux_l[1])/dx + (flux_u[1] - flux_d[1])/dr)*DT;
                            flux[2] = ((flux_r[2] - flux_l[2])/dx + (flux_u[2] - flux_d[2])/dr)*DT;
                            flux[3] = ((flux_r[3] - flux_l[3])/dx + (flux_u[3] - flux_d[3])/dr)*DT;
                        // }

                        /* SOURCE TERM CORRESPONDING TO CYLINDRICAL SYSTEM OF COORDINATES */ 
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

                                flux[0] -= source[0]*DT;
                                flux[1] -= source[1]*DT;
                                flux[2] -= source[2]*DT;
                                flux[3] -= source[3]*DT;
                                // printf("%f %f %f %f \n",source[0],source[1],source[2],source[3]);

                                if (isnan(source[0]) || isnan(source[1]) || isnan(source[2]) || isnan(source[3]))
                                {
                                    printf("DOM UP IS %d\n", dom_up);
                                    fprintf(stderr,"\nERROR:SOURCE IS NaN FOR VALVE %d at time %d and at domain %d %d %d: %f %f %f %f (ytip is %f)!\n",m,t,k,i,j,source[0],source[1],source[2],source[3],y_tip);
                                    exit(0);
                                }
                            }
                        }

                        /* INCREMENT RUNGE-KUTTA VARIABLES WITH CORRESPONDING FLUX */
                        RK_k = 1./(RK_ORDER-trk);
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

        append_data_to_file(p_wall,DT*t,P_MEAN_WALL_FILENAME,EXP_EXTENSION);
        append_data_to_file(mfr_intake,DT*t,INTAKE_MFR_FILENAME,EXP_EXTENSION);
        append_data_to_file(p_tube,DT*t,P_TUBE_FILENAME,EXP_EXTENSION);
        append_data_to_file(rho_tube,DT*t,RHO_TUBE_FILENAME,EXP_EXTENSION);

        if (PLENUM_ON==1)
        {
            mfr_plenum = 0.0, p_plenum = 0.0, rho_plenum = 0.0, p_drag = 0.0;

            mfr_plenum = mfr_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],u[dom_up],NGHOST,NGHOST);
            p_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NGHOST);
            rho_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],NGHOST);
            p_drag = average_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NXtot[dom_up]-NGHOST-1,NGHOST);

            append_data_to_file(p_plenum,DT*t,PLEN_P_FILENAME,EXP_EXTENSION);
            append_data_to_file(rho_plenum,DT*t,PLEN_RHO_FILENAME,EXP_EXTENSION);
            append_data_to_file(mfr_plenum,DT*t,PLENUM_MFR_FILENAME,EXP_EXTENSION);
            append_data_to_file(p_drag,DT*t,PLENUM_DRAG_PRESSURE_FILENAME,EXP_EXTENSION);
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
            append_data_to_file(total_mfr,DT*t,MFR_TOT_FILENAME,EXP_EXTENSION);
            append_data_to_file(total_pfr,DT*t,PFR_FILENAME,EXP_EXTENSION);
            append_multidata_to_file(N_VALVE,stage_mfr,DT*t,MFR_FILENAME,EXP_EXTENSION);

            // Valve displacement and pressure ratio export
            append_multidata_to_file(N_VALVE,ytip,DT*t,Y_TIP_FILENAME,EXP_EXTENSION);
            append_multidata_to_file(N_VALVE,pratio,DT*t,P_RATIO_FILENAME,EXP_EXTENSION);
        }

        /* EXPORT ALL FLUID AND SOLID DATA EVERY N STEPS */
        if(t%N_EXPORT==0)
        {
            // Fluid data export
            export_fluid_data(NDOMAIN,NXtot,NYtot,x,y,xc,yc,rho,u,v,p,E,T,H,DT*t,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT);

            // Valve data export
            if (SOLID_ON==1)
            {
                // Valve displacement
            	export_valve_data(N_VALVE,N_FEM+1,x_FEM,y_FEM,R0,DT*t,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT,NGHOST);
            }
        }


        /* DISPLAY CFL */
        if(t%N_CFL==0)
        {
            // COMPUTE CFL (BASED ON TOTAL VELOCITY)
            CFL = get_cfl(x,y,u,v,T,R,GAMMA,NDOMAIN,NXtot,NYtot,DT);

            // COUNTING TIME BETWEEN ITERATIONS
            clock_gettime(CLOCK_REALTIME,&end_time_loop);
            SIM_TIME_LOOP = (end_time_loop.tv_sec - start_time_loop.tv_sec)+(end_time_loop.tv_nsec - start_time_loop.tv_nsec)/1E9;
            printf("Iteration %d (%.1f %%)... Max. CFL at t=%.5f sec is: %f. \n[completed in %f sec]\n",t,(double)(t)/Ntstep*100.0,DT*t,CFL,SIM_TIME_LOOP);
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

