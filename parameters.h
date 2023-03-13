
/* PARAMETERS TO BE SPECIFIED BY THE USER TO RUN A SIMULATION */
/*
This file needs to be adapted to the case the user wants to run (number of domains, use of previously obtained numerical data, choice of initial solution, etc.) so be careful what you do!
*/


/* SIMULATION CASES DEFINED BY PREVIOUS USERS */
// char *SIM_CASE = "det_tube"; // A Sod shock tube case solved in 2D
// char *SIM_CASE = "mw_tube"; // A simple MW rocket and intake region without reed valve
char *SIM_CASE = "plenum_rocket"; // A simple MW rocket, intake, and a plenum region
// char *SIM_CASE = "ground_rocket"; // A simple MW rocket and intake
// char *SIM_CASE = "sup_plenum"; // A single plenum region in a highly supersonic flow
// char *SIM_CASE = "valve_tank"; // Two tanks connected by a valve
// char *SIM_CASE = "sup_plen_rocket"; // A single plenum region in a highly supersonic flow


/* PRE-RUN CASE (IF NEEDED) */
// char *PRERUN_FOLDER = "PRERUN/M1.968_PL0.0070_osc/";
// char *PRERUN_FOLDER = "PRERUN/M1.968_PL0.0090_osc/";
// char *PRERUN_FOLDER = "PRERUN/M1.968_PL0.0200_osc/";
// char *PRERUN_FOLDER = "PRERUN/M1.968_PL0.0270_osc/";
// char *PRERUN_FOLDER = "PRERUN/M1.968-0.0184/";
// char *PRERUN_FOLDER = "PRERUN/M1.968_PL0.0110/";
// char *PRERUN_FOLDER = "PRERUN/M1.968-PL0.0160/";
char *PRERUN_FOLDER = "PRERUN/last/";


/* GENERAL FLUID PARAMETERS */ 
double R           = 287.0;    // R constant for gas state equation
double GAMMA       = 1.4;      // Heat capacity ratio (constant)


/* GEOMETRY PARAMETERS */
double L0          = 0.01;      // Length of reference domain 0
double R0          = 0.028;//0.125;//0.028;     // Height of reference domain 0
double L_OUT       = 4.0;//8.0;//4.0;       // Outlet length    
double R_OUT       = 4.0;//8.0;//4.0;       // Outlet above-tube height to ref. height
double L_TUBE      = 0.5;//1.0;//0.397;     // Length of the tube
double H_WALL      = 0.002;     // Wall thickness of the thruster


/* VALVE GEOMETRY */
double B0         = 16.0E-3;   // Width of valve at fixed end
double B1         = 12.0E-3;   // Width of valve at free end
double H0         = 0.35E-3;   // Thickness of valve at fixed end
double H1         = 0.15E-3;   // Thickness  of valve at free end
double L_V        = 25.0E-3;   // Length of valve ( and of corresponding fluid domain)
double L_FIX      = 6.0E-3;    // Additional length of reed plate that is used for fixation 
double E_V        = 110.0E9;   // Young modulus (of TiAl alloy)
double RHO_V      = 4400.0;    // Density of valve
double F0         = 470.0;     // Natural frequency (Hz) of reed petal (approximate)
int N_FEM         = 30;        // Number of FEM elements
int N_FIX         = 3;         // Number of nodes over L_FIX length
int N_CLAMP       = 2;         // Number of clamped nodes
int N_V_PER_STAGE = 8;         // Number of reed valves per stage
double L_HOLE     = 23.0E-3;   // Length of the hole through which flow runs from plenum to thruster
double B_HOLE     = 10.0E-3;   // Width of hole through which flow runs from plenum to thruster
double HOLE_FACTOR = 0.5;      // Proportion of hole length over which mass flow rate is injected


/* PLENUM GEOMETRY */
double L_P       = 0.498;                           // Length of plenum region (remove wall length)
double H_P       = 0.011;                           // Height of plenum region (above thruster exterior wall)
double H_P_OUT   = 0.55;//0.45;                           // Free-stream region above plenum's height
double L_INLET_P = 3.5;//0.10;//3.5;//0.20;//3.5;         // Inlet distance with regard to plenum most upstream point
double NOSE_ANGLE = 13.98*M_PI/180.0;               // Angle in radian of cone angle for the rocket

/* FINITE ELEMENT MODEL */
int N_DOF_PER_NODE    = 2;         // Number of DOFs per node (beam element in 2D)
double RAYLEIGH_ALPHA = 0.0;       // Alpha coef. for Rayleigh damping (alpha*M + beta*K)
double RAYLEIGH_BETA  = 5.0e-6;    // Beta coef. for Rayleigh damping
double C1 = 5.0E-8;//5.0E-7        // C1 damping factor in: eps_eff = C1 + y_tip_dot*C2 (y_tip_dot>0)
double C2 = 0.0E-8;//2.0E-8
double C3 = 0.0007;                // Flow damping coefficient for closing: eps_eff = C3*y_tip (y_tip_dot<0)
// double YTIP_RATIO = 2.0;


/* TIME-RELATED PARAMETERS */
double TSIM        = 0.01;        // Simulation time in second (Final time)
double DT          = 10.0E-7;      // Time step in second
char *W_FORMAT     = "%.10g";     // Adapt writing format to time step value
int N_STEP_CFL     = 200;         // Maximum number of displayed CFL time steps
int N_STEP_EXP     = 200;         // Maximum number of exported time steps
int N_STEP_P_MFR   = 500;         // Maximum number of exported pressure/MFR steps


/* CARTESIAN GRID FOR DEFAULT REGION (REGION 0) */
int NX                 = 250;//250;    // Number of cells along X over domain 0
int NY                 = 14;//30;     // Number of cells along Y over domain 0
double MESH_RATIO_OUT  = 8.0;    // Mesh size ratio between above/below valve regions
double MESH_RATIO_IN   = 8.0;    // Mesh size ratio in the above plenum region


/* INITIAL AMBIENT CONDITIONS & MICROWAVE BEAMING */
double P0          = 101325.0;   // Ambient pressure (Pa)
double T0          = 298.0;      // Ambient temperature (K)
double S0          = 2000.0E3;   // Total microwave beam power (W)
double ETA         = 1.0;        // Energy absorption coefficient


/* ALGORITHM PARAMETERS (CHANGE TO YOUR OWN RISK) */
int NGHOST         = 2;      // Number of ghost cells added at boundaries (MUSCL)
double MUSCL_BIAS  = 1.0/3;  // MUSCL bias coefficient between upwind and downwind differences
int LIMITERNAME    = 0;      // Flux Limiter function for MUSCL (-1: none / 0:minmod / 1:superbee / 2:vanAlbada1 / 3:vanAlbada2)
double AUSM_K      = 10.0;   // Bias parameter for AUSM switching function
double ENTRO_FIX_C = 0.125;  // Parameter for entropy fix in AUSM-DV scheme
int RK_ORDER       = 4;      // Order of Runge-Kutta time-iteration scheme


/* EXPORT FORMAT, FOLDER AND PARAMETERS */
char *EXP_EXTENSION   = ".dat";                        // Format of exported files
char *OUT_FOLDERNAME = "OUTPUT";                      // Name of directory in which results are exported
char *PAR_FILENAME   = "OUTPUT/Parameters.dat";       // Filename of the parameter file

/* USER-DEFINED EXPORT VARIABLES */ 
char *P_TUBE_FILENAME      = "OUTPUT/P_tube_mean";       // Average pressure in tube
char *RHO_TUBE_FILENAME    = "OUTPUT/RHO_tube_mean";     // Average density in tube
char *P_MEAN_WALL_FILENAME = "OUTPUT/p_wall";            // Mean wall pressure history at wall (cyl.)
char *PLENUM_DRAG_PRESSURE_FILENAME = "OUTPUT/p_drag";   // Mass flow rate at plenum intake
char *Y_TIP_FILENAME       = "OUTPUT/y_tip";             // Valve tip displacement history
char *MFR_TOT_FILENAME     = "OUTPUT/mfr_total";         // Total mass flow rate history through valves
char *MFR_FILENAME         = "OUTPUT/mfr_stage";         // Mass flow rate at each reed stage
char *INTAKE_MFR_FILENAME  = "OUTPUT/mfr_intake";        // Intake mass flow rate
char *PLENUM_MFR_FILENAME  = "OUTPUT/mfr_plenum";        // Mass flow rate at plenum intake
char *P_RATIO_FILENAME     = "OUTPUT/p_ratio";           // Pressure ratio history at each valve
char *PFR_FILENAME         = "OUTPUT/pfr";               // Partial Filling Rate history
char *PLEN_P_FILENAME      = "OUTPUT/plenum_pressure";   // Intake mass flow rate
char *PLEN_RHO_FILENAME    = "OUTPUT/plenum_density";    // Intake mass flow rate 

/* FEM MATRICES FILENAMES */
char *FEM_SECTIONS = "OUTPUT/FEM";      // Sectional information of the FEM model
char *FEM_DOF      = "OUTPUT/DOF";      // List of degrees of freedom
char *FEM_K        = "OUTPUT/K";        // FEM Stiffness matrix
char *FEM_M        = "OUTPUT/M";        // FEM Mass matrix
char *FEM_C        = "OUTPUT/C";        // FEM Reynolds Damping matrix

/* Attribution of the values inside the boundary condition array. */
char *** init_boundary_type(int Ndom, char *boundary_type[][4])
{
    int char_len = 4;
    char ***bounds = malloc(Ndom*4*sizeof(char)*char_len);
    for (int k = 0; k < Ndom; ++k)
    {
        bounds[k] = malloc(char_len*4*sizeof(char));
        for (int j = 0; j < 4; ++j)
        {
            bounds[k][j] = malloc(char_len*sizeof(char));
            strcpy(bounds[k][j],boundary_type[k][j]);
        }
    }
    return bounds;
}
/* Memory allocation of the boundary condition array. */
char * init_boundary_location(char boundary_location[4])
{
    char *bounds = malloc(4*sizeof(char));
    for (int i = 0; i < 4; ++i)
    {
        bounds[i] = boundary_location[i];
    }

    return bounds;
}

/* Initialization procedure of parameters depending on case name specified by the user. New IF loops should be added if the user decides to add new cases. */
void init_case(char *case_name, double *l_v_tot, double *M_ref, int *lower_domain, int *upper_domain,double **x0, double **y0, double **length_x, double **length_y, double **x_ratio, double **y_ratio, double **xv0, char **bound_loc, char ****bound_type, int **nx, int **ny, int *ndom, int *nv, int *solid_on, int *plenum_on)
{
    /* Common to all cases: boundary definition */
    char BOUNDARY_ORDER[]  = {'L'   ,'R'   ,'D'   ,'U'};

    if (strcmp(case_name,"det_tube")==0)
    {
        *ndom = 2;
        *nv = 0;
        *solid_on = 0;
        *l_v_tot = L_V + L_FIX;
        *M_ref = 0.0;
        *lower_domain = 0;
        *upper_domain = 0;

        double xstart[] = {0.0,L_TUBE}; // X location of most bottom left point      
        double ystart[] = {0.0,0.0}; // Y location of most bottom left point
        double xlength[] = {L_TUBE,L_TUBE}; // X-wise length of each domain
        double ylength[] = {R0,R0}; // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},  // Boundaries of domain 0 (Tube)
                                    {"con0","slip","slip","slip"}}; 

        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,1.0};
        double graded_ratio_y[] = {1.0,1.0};

        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }
    }
    else if (strcmp(case_name,"mw_tube")==0)
    {
        *ndom = 6;
        *nv = 0;
        *solid_on = 0;
        *l_v_tot = L_V + L_FIX;
        *M_ref = 0.0;
        *lower_domain = 0;
        *upper_domain = 4;

        // double xstart[] = {0.0,L_TUBE,L_TUBE}; // X location of most bottom left point      
        // double ystart[] = {0.0,0.0,R0}; // Y location of most bottom left point
        // double xlength[] = {L_TUBE,L_OUT,L_OUT}; // X-wise length of each domain
        // double ylength[] = {R0,R0,R_OUT}; // Y-wise length of each
        // char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},  // Boundaries of domain 0 (Tube)
        //                             {"con0","slip","slip","con2"},  // Boundaries of domain 2 (Outlet 1)
        //                             {"slip","slip","con1","slip"}};  // Boundaries of domain 2 (Outlet 2)
        // double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        // double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        // double graded_ratio_x[] = {1.0,MESH_RATIO_OUT,MESH_RATIO_OUT};
        // double graded_ratio_y[] = {1.0,1.0,MESH_RATIO_OUT};

        double xstart[] = {0.0,L_TUBE,L_TUBE,-L_INLET_P,0.0,L_TUBE}; // X location of most bottom left point      
        double ystart[] = {0.0,0.0,R0,R0+H_WALL,R0+H_WALL,R0+H_WALL}; // Y location of most bottom left point
        double xlength[] = {L_TUBE,L_OUT,L_OUT,L_INLET_P,L_TUBE,L_OUT}; // X-wise length of each domain
        double ylength[] = {R0,R0,H_WALL,R_OUT,R_OUT,R_OUT}; // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},  // Boundaries of domain 0 (Tube)
                                    {"con0","slip","slip","con2"},  // Boundaries of domain 2 (Outlet 1)
                                    {"slip","slip","con1","con5"},  // Boundaries of domain 2 (Outlet 2)
                                    {"slip","con4","slip","slip"},  // Boundaries of domain 3 (Plenum inlet)
                                    {"con3","con5","slip","slip"},  // Boundaries of domain 4 (Plenum)
                                    {"con4","slip","con2","slip"}}; // Boundaries of domain 4 (Facing plenum wall)

        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_IN,1.0,MESH_RATIO_OUT};
        double graded_ratio_y[] = {1.0,1.0,1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_OUT};

        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }
    }
    else if (strcmp(case_name,"plenum_rocket")==0)
    {
        *ndom = 8;
        *nv = 6;
        *solid_on = 1;
        *l_v_tot = L_V + L_FIX;
        *M_ref = 0.0;
        *lower_domain = 0;
        *upper_domain = 4;

        double xstart[] = {0.0,L_TUBE,L_TUBE,-L_INLET_P,0.0,-L_INLET_P,-L_INLET_P,0.0};// X location of most bottom left point      
        double ystart[] = {0.0,0.0,R0,R0+H_WALL,R0+H_WALL,R0+H_WALL+H_P,R0+3.0*H_WALL+H_P,R0+3.0*H_WALL+H_P};        // Y location of most bottom left point
        double xlength[] = {L_TUBE,L_OUT,L_OUT,L_INLET_P,L_P,L_INLET_P,L_INLET_P,L_TUBE-H_WALL}; // X-wise length of each domain
        double ylength[] = {R0,R0,R_OUT,H_P,H_P,2.0*H_WALL,H_P_OUT,H_P_OUT};          // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},  // Boundaries of domain 0 (Tube)
                                    {"con0","slip","slip","con2"},  // Boundaries of domain 2 (Outlet 1)
                                    {"slip","slip","con1","slip"},  // Boundaries of domain 2 (Outlet 2)
									{"slip","con4","slip","con5"},  // Boundaries of domain 3 (Plenum inlet)
                                    {"con3","slip","slip","slip"},  // Boundaries of domain 4 (Plenum)
                                    {"slip","slip","con3","con6"},  // Boundaries of domain 4 (Facing plenum wall)
                                    {"slip","con7","con5","slip"},  // Boundaries of domain 5 (Above plenum inlet)
                                    {"con6","slip","slip","slip"}}; // Boundaries of domain 6 (Above plenum)
        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_OUT,1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,1.0};
        double graded_ratio_y[] = {1.0,1.0,MESH_RATIO_OUT,1.0,1.0,1.0,MESH_RATIO_IN,MESH_RATIO_IN};
        double x_v_start[] = {0.04,0.081,0.129,0.170,0.218,0.259};

        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }

        double *X_V_START = malloc(*nv*sizeof(double));
        *xv0 = X_V_START;
        for (int i = 0; i < *nv; ++i)
        {
        	X_V_START[i] = x_v_start[i];
        }
    }
    else if (strcmp(case_name,"ground_rocket")==0)
    {
        *ndom = 6;
        *nv = 6;
        *solid_on = 1;
        *l_v_tot = L_V + L_FIX;
        *M_ref = 0.0;
        *lower_domain = 0;
        *upper_domain = 4;

        double xstart[] = {0.0,L_TUBE,L_TUBE,-L_INLET_P,0.0,L_TUBE}; // X location of most bottom left point      
        double ystart[] = {0.0,0.0,R0,R0+H_WALL,R0+H_WALL,R0+H_WALL}; // Y location of most bottom left point
        double xlength[] = {L_TUBE,L_OUT,L_OUT,L_INLET_P,L_TUBE,L_OUT}; // X-wise length of each domain
        double ylength[] = {R0,R0,H_WALL,R_OUT,R_OUT,R_OUT}; // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},  // Boundaries of domain 0 (Tube)
                                    {"con0","slip","slip","con2"},  // Boundaries of domain 2 (Outlet 1)
                                    {"slip","slip","con1","con5"},  // Boundaries of domain 2 (Outlet 2)
                                    {"slip","con4","slip","slip"},  // Boundaries of domain 3 (Plenum inlet)
                                    {"con3","con5","slip","slip"},  // Boundaries of domain 4 (Plenum)
                                    {"con4","slip","con2","slip"}}; // Boundaries of domain 4 (Facing plenum wall)
        // *ndom = 5;
        // *nv = 6;
        // *solid_on = 1;
        // *l_v_tot = L_V + L_FIX;
        // *M_ref = 0.0;
        // *lower_domain = 0;
        // *upper_domain = 4;

        // double xstart[] = {0.0,L_TUBE,L_TUBE,-L_INLET_P,0.0}; // X location of most bottom left point      
        // double ystart[] = {0.0,0.0,R0,R0+H_WALL,R0+H_WALL}; // Y location of most bottom left point
        // double xlength[] = {L_TUBE,L_OUT,L_OUT,L_INLET_P,L_TUBE-H_WALL}; // X-wise length of each domain
        // double ylength[] = {R0,R0,R_OUT,R_OUT,R_OUT}; // Y-wise length of each
        // char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},  // Boundaries of domain 0 (Tube)
        //                             {"con0","slip","slip","con2"},  // Boundaries of domain 2 (Outlet 1)
        //                             {"slip","slip","con1","slip"},  // Boundaries of domain 2 (Outlet 2)
        //                             {"slip","con4","slip","slip"},  // Boundaries of domain 3 (Plenum inlet)
        //                             {"con3","slip","slip","slip"}};  // Boundaries of domain 4 (Plenum)

        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_IN,1.0,MESH_RATIO_OUT};
        double graded_ratio_y[] = {1.0,1.0,1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_OUT};
        // double graded_ratio_x[] = {1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_IN,1.0};
        // double graded_ratio_y[] = {1.0,1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,MESH_RATIO_OUT};
        double x_v_start[] = {0.04,0.079,0.129,0.168,0.218,0.257,0.307,0.346}; // exp
        // double x_v_start[] = {0.04,0.2,0.129,0.168,0.218,0.257,0.307,0.346}; // gradually added
        // double x_v_start[] = {0.001,0.067,0.133,0.199,0.265,0.331}; // distrib
        // double x_v_start[] = {0.001,0.038,0.075,0.112,0.149,0.186}; // wall
        // double x_v_start[] = {0.001,0.041,0.081,0.121,0.161,0.201,0.241,0.281}; // customizable
        // double x_v_start[] = {L_TUBE*0.5};

        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }

        double *X_V_START = malloc(*nv*sizeof(double));
        *xv0 = X_V_START;
        for (int i = 0; i < *nv; ++i)
        {
            X_V_START[i] = x_v_start[i];
        }
    }
    else if (strcmp(case_name,"sup_plenum")==0)
    {
        *ndom = 5;
        *solid_on = 0;
        *lower_domain = 0;
        *upper_domain = 1;
        *plenum_on = 1;

        // 40 km and 2 km/s
        // P0 = 277.546;
        // T0 = 251;
        // *M_ref = 6.30;
        
        // // 10 km and 0.6 km/s
        // P0 = 26437;
        // T0 = 223;
        // *M_ref = 2.0;

        // 10 km and 0.6 km/s, shock angle of 13.98 degrees and isentropic expansion fan
        P0 = 26569.0;
        T0 = 226.0;
        *M_ref = 1.968;

        // 10 km and 1 km/s
        // P0 = 26437;
        // T0 = 223;
        // *M_ref = 3.34;

        *l_v_tot = L_V + L_FIX;

        // double l_nose = (R0+H_WALL)/tan(NOSE_ANGLE);
        // *cone_index = 1;
        // *precone_index = 0;

        // // X location of most bottom left point
        // double xstart[] = { -L_INLET_P , -l_nose , 0.0 , -L_INLET_P , -l_nose , 0.0 };

        // // Y location of most bottom left point
        // double ystart[] = {  R0+H_WALL , R0+H_WALL , R0+H_WALL , R0+H_WALL+H_P , R0+H_WALL+H_P , R0+H_WALL+H_P };
        // double xlength[] = { L_INLET_P - l_nose , l_nose , L_P , L_INLET_P - l_nose , l_nose , L_P }; // X-wise length of each domain
        // double ylength[] = { H_P , H_P , H_P , H_P_OUT , H_P_OUT , H_P_OUT };          // Y-wise length of each
        // char *BOUNDARY_TYPE[][4] = {{"supI","con1","slip","con3"},  // Boundaries of domain 0
        //                             {"con0","con2","incl","con4"},  // Boundaries of domain 2
        //                             {"con1","slip","slip","slip"},  // Boundaries of domain 2
        //                             {"supI","con4","con0","slip"},
        //                             {"con3","con5","con1","slip"},
        //                             {"con4","supO","slip","slip"},
                                    // }; // Boundaries of domain 2

                // X location of most bottom left point
        double xstart[] = { -L_INLET_P, 0.0 , -L_INLET_P, -L_INLET_P , 0.0};

        // Y location of most bottom left point
        double ystart[] = {  R0+H_WALL , R0+H_WALL , R0+H_WALL+H_P, R0+2.0*H_WALL+H_P , R0+2.0*H_WALL+H_P};
        double xlength[] = { L_INLET_P, L_P , L_INLET_P, L_INLET_P , L_P}; // X-wise length of each domain
        double ylength[] = { H_P , H_P , H_WALL , H_P_OUT , H_P_OUT};          // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"supI","con1","slip","con2"},  // Boundaries of domain 0
                                    {"con0","slip","slip","slip"},
                                    {"supI","slip","con0","con3"},
                                    {"supI","con4","con2","supO"},
                                    {"con3","supO","slip","supO"}};

        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,1.0,1.0,1.0,1.0};
        double graded_ratio_y[] = {1.0,1.0,1.0,1.0,1.0};

        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }
    }
    else if (strcmp(case_name,"sup_plen_rocket")==0)
    {
        *ndom = 8;
        *nv = 6;
        *solid_on = 1;
        *lower_domain = 0;
        *upper_domain = 4;
        *plenum_on = 1;

        // 40 km and 2 km/s
        // P0 = 277.546;
        // T0 = 251;
        // *M_ref = 6.30;
        
        // 10 km and 0.5 km/s
        P0 = 26569.0;
        T0 = 226.0;
        *M_ref = 1.968;


        *l_v_tot = L_V + L_FIX;

        double xstart[] = {0.0,L_TUBE,L_TUBE,-L_INLET_P,0.0,-L_INLET_P,-L_INLET_P,0.0};         // X location of most bottom left point      
        double ystart[] = {0.0,0.0,R0,R0+H_WALL,R0+H_WALL,R0+H_WALL+H_P,R0+2.0*H_WALL+H_P,R0+2.0*H_WALL+H_P};        // Y location of most bottom left point
        double xlength[] = {L_TUBE,L_OUT,L_OUT,L_INLET_P,L_P,L_INLET_P,L_INLET_P,L_P}; // X-wise length of each domain
        double ylength[] = {R0,R0,R_OUT,H_P,H_P,H_WALL,H_P_OUT,H_P_OUT};          // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"slip","con1","slip","slip"},
                                    {"con0","slip","slip","con2"},
                                    {"slip","slip","con1","slip"},
                                    {"supI","con4","slip","con5"},
                                    {"con3","slip","slip","slip"},
                                    {"supI","slip","con3","con6"},
                                    {"supI","con7","con5","slip"},
                                    {"con6","supO","slip","slip"}};
        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,MESH_RATIO_OUT,MESH_RATIO_OUT,1.0,1.0,1.0,1.0,1.0};
        double graded_ratio_y[] = {1.0,1.0,MESH_RATIO_OUT,1.0,1.0,1.0,1.0,1.0};
        // double x_v_start[] = {0.04,0.081,0.129,0.170,0.218,0.259};
        double x_v_start[] = {0.04,0.079,0.129,0.168,0.218,0.257}; // exp
        // double x_v_start[] = {0.04,0.081,0.129,0.170,0.218,0.259,0.307,0.348,0.396,0.437};


        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }

        double *X_V_START = malloc(*nv*sizeof(double));
        *xv0 = X_V_START;
        for (int i = 0; i < *nv; ++i)
        {
            X_V_START[i] = x_v_start[i];
        }
    }
    else if (strcmp(case_name,"valve_tank")==0)
    {
        *ndom = 2;
        // *nv = 1;
        *nv = 1;
        *solid_on = 1;
        *l_v_tot = L_V + L_FIX;
        *M_ref = 0.0;
        *lower_domain = 0;
        *upper_domain = 1;

        double xstart[] = {0.0,0.0};        // X location of most bottom left point      
        double ystart[] = {0.0,R0+H_WALL};  // Y location of most bottom left point
        double xlength[] = {L_TUBE,L_TUBE}; // X-wise length of each domain
        double ylength[] = {R0,R0};          // Y-wise length of each
        char *BOUNDARY_TYPE[][4] = {{"slip","slip","slip","slip"},  // Boundaries of domain 0 (Tube)
                                    {"slip","slip","slip","slip"}}; // Boundaries of domain 2 (Outlet 1)
        double dx = L_TUBE/NX; // Basic X-wise grid size used to estimate number of cells
        double dy = R0/NY;  // Basic Y-wise grid size used to estimate number of cells
        double graded_ratio_x[] = {1.0,1.0};
        double graded_ratio_y[] = {1.0,1.0};
        // double x_v_start[] = {0.15};
        double x_v_start[] = {0.18,0.08,0.13,0.18,0.23,0.28};

        double *XSTART = malloc(*ndom*sizeof(double));
        double *YSTART = malloc(*ndom*sizeof(double));
        double *XLENGTH = malloc(*ndom*sizeof(double));
        double *YLENGTH = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_X = malloc(*ndom*sizeof(double));
        double *GRID_RATIO_Y = malloc(*ndom*sizeof(double));
        int *NXtot = malloc(*ndom*sizeof(int));
        int *NYtot = malloc(*ndom*sizeof(int));

        char *B_LOC = init_boundary_location(BOUNDARY_ORDER);
        char ***B_TYPE = init_boundary_type(*ndom,BOUNDARY_TYPE);

        *x0 = XSTART;
        *y0 = YSTART;
        *length_x = XLENGTH;
        *length_y = YLENGTH;
        *x_ratio = GRID_RATIO_X;
        *y_ratio = GRID_RATIO_Y;
        *nx = NXtot;
        *ny = NYtot;
        *bound_loc = B_LOC;
        *bound_type = B_TYPE;

        for (int i = 0; i < *ndom; ++i)
        {
            XSTART[i] = xstart[i];
            YSTART[i] = ystart[i];
            XLENGTH[i] = xlength[i];
            YLENGTH[i] = ylength[i];

            GRID_RATIO_X[i] = graded_ratio_x[i];
            GRID_RATIO_Y[i] = graded_ratio_y[i];

            NXtot[i] = (int)fmax(2+2*NGHOST,(XLENGTH[i]/dx/graded_ratio_x[i]+2*NGHOST));
            NYtot[i] = (int)fmax(2+2*NGHOST,(YLENGTH[i]/dy/graded_ratio_y[i]+2*NGHOST));
        }

        double *X_V_START = malloc(*nv*sizeof(double));
        *xv0 = X_V_START;
        for (int i = 0; i < *nv; ++i)
        {
            X_V_START[i] = x_v_start[i];
        }
    }
    else
    {
        fprintf(stderr, "\nERROR: Specified case was not recognized! \n"); exit(0);
    }
}

/* End of file */