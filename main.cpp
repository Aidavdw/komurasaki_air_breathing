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
//#include "case_det_tube.cpp" // Temporary hard-code of specific case implementation
#include <chrono>

#include "reed_valve.h"


/* MAIN ROUTINE */
int main()
{
    /* STARTING PROGRAM, GETTING TIME... */
    std::cout << "Program starting...\n";
    std::time_t timeAtStartOfProgram = std::time(0);

    // INITIALIZE CASE
    SimCase simCase;
    //LoadCase(&simCase);

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

    // Todo: make actual thing that handles how reed valves etc are chosen
    ReedValve reedValve = ReedValve();
    simCase.RegisterValve(reedValve);


    /* INITIAL CONDITIONS ON DOMAINS */
    simCase.ApplyInitialConditions();
    std::cout << "Initial conditions applied." << std::endl;






    /* TIME LOOP VARIABLES */
    double CFL=0;                        // For display of CFL
    double SIM_TIME_LOOP;                // To monitor real time between two exports
    double total_mfr = 0.0, total_pfr = 0.0;


    /* COUNTING TIME BETWEEN ITERATIONS */
    auto timeAtStartOfCalculations = std::chrono::system_clock::now();

    std::cout << "Starting time loop." << std::endl;
    /* START TIME LOOP */
    for (int timeStepNumber = 1; timeStepNumber < simCase.totalSimulationTimeStepCount; ++timeStepNumber)
    {
        /* 4TH ORDER RUNGE-KUTTA PREDICTOR-CORRECTOR LOOP (iterate 4 times)*/
        for (int rungeKuttaIterationNumber = 0; rungeKuttaIterationNumber < RK_ORDER; ++rungeKuttaIterationNumber)
        {
            /* FOR FIRST ITERATION, TAKE PARAMETERS OF PREVIOUS SOLUTION */
            for (auto& domainIter : simCase.domains)
            {
                Domain& domain = domainIter.second;

                //#pragma omp parallel for // START PARALLEL LOOP
                for (int xIndex = 0; xIndex < domain.size[0]; ++xIndex)
                {
                    for (int yIndex = 0; yIndex < domain.size[1]; ++yIndex)
                    {
                        // Take solution from previous iteration for this new iteration
                        if(rungeKuttaIterationNumber==0)
                        {
                            // First Runge-Kutta iteration
                            domain.CopyFieldQuantitiesToBuffer(EFieldQuantityBuffer::MAIN, EFieldQuantityBuffer::RUNGEKUTTA);
                        }
                        else
                        {
                            // Second, third, etc. iterations...
                            domain.CopyFieldQuantitiesToBuffer(EFieldQuantityBuffer::T, EFieldQuantityBuffer::RUNGEKUTTA);
                        }

                        // Reset all variable arrays used for final memory saving
                        // A: if they're being written further down, no need to set them to 0, they might just be overwritten.
                        rhot[domainNumber][xIndex][yIndex]=0;
                        ut[domainNumber][xIndex][yIndex]=0;
                        vt[domainNumber][xIndex][yIndex]=0;
                        Et[domainNumber][xIndex][yIndex]=0;
                        pt[domainNumber][xIndex][yIndex]=0;
                        Tt[domainNumber][xIndex][yIndex]=0;
                        Ht[domainNumber][xIndex][yIndex]=0;

                        // In addition, reset "sonicPoints" arrays
                        sonic_x[domainNumber][xIndex][yIndex]=0;
                        sonic_y[domainNumber][xIndex][yIndex]=0;
                    }
                }
            }

            for (IValve& valve : simCase.valves)
            {
                valve.Update();
            }
            
            //#pragma omp parallel for num_threads(N_VALVE)
            for (int valveIndex = 0; valveIndex < N_VALVE; ++valveIndex)
            {
                // Compute pressure difference at FEM nodes and deduce load on each element/node
                // fem_pressure(fem_n[domainNumber],fem_index_inf[domainNumber],x_FEM[domainNumber],NYtot[dom_low],n_cell_p_fem,pRK[dom_low],pRK[dom_up],NGHOST,dp_interface[domainNumber]);
                // fem_load(N_FEM,N_DOF,N_DOF_PER_NODE,N_CLAMP,NGHOST,X_V_START[domainNumber],p_neighbour[domainNumber],fem_n[domainNumber],dp_interface[domainNumber],x[dom_low],x_FEM[domainNumber],B0,B1,L_T,F_DOF[domainNumber]);

                fem_pressure(N_NODE,fem_n[valveIndex],fem_index_inf[valveIndex],x_FEM[valveIndex],p_neighbour[valveIndex],p_coef[valveIndex],NYtot[dom_low],pRK[dom_low],pRK[dom_up],NGHOST,p_FEM[valveIndex]);
                fem_load(N_FEM,N_DOF,N_CLAMP,N_DOF_PER_NODE,p_FEM[valveIndex],U2_DOF[valveIndex],F_DOF[valveIndex],x_FEM[valveIndex],b);
                fem_flow_damping(N_FEM,N_DOF,N_CLAMP,N_DOF_PER_NODE,U1_DOF[valveIndex],U2_DOF[valveIndex],b,h,RHO_V,F0,DT,C1,C2,C3,F_DOF[valveIndex]);

                // Solve FEM system and obtain Runge-Kutta valve distorsion at current RK loop
                newmark_solve(N_DOF,N_ACTIVE,L_R1,R2,R3,F_DOF[valveIndex],act_DOF,U1_DOF[valveIndex],U2_DOF[valveIndex],U2_DOF_K[valveIndex]);

                // Update FEM mesh and contrain displacement in positive domain
                update_valve(N_NODE,N_DOF_PER_NODE,U0_DOF[valveIndex],U1_DOF[valveIndex],U2_DOF_K[valveIndex],y_FEM[valveIndex]);
            }

            if (SOLID_ON==1)
            {
                // #pragma omp parallel for num_threads(N_VALVE) (bug if activated?)
                for (int valveIndex = 0; valveIndex < N_VALVE; ++valveIndex)
                {
                    // Compute mean pressure on each side and mass-flow rate at valve
                    mean_p_inf[valveIndex] = mean_at_valve(fem_index_inf[valveIndex],fem_n[valveIndex],NYtot[dom_low]-NGHOST-n_cell_p_fem,n_cell_p_fem,pRK[dom_low]);
                    mean_p_sup[valveIndex] = mean_at_valve(fem_index_inf[valveIndex],fem_n[valveIndex],NGHOST,n_cell_p_fem,pRK[dom_up]);
                    mean_rho_sup[valveIndex] = mean_at_valve(fem_index_inf[valveIndex],fem_n[valveIndex],NGHOST,n_cell_p_fem,rhoRK[dom_up]);

                    // Tip deflection and angle
                    double y_tip = y_FEM[valveIndex][N_FEM-1];

                    if (isnan(mean_p_sup[valveIndex]) || isnan(mean_p_inf[valveIndex]) || isnan(mean_rho_sup[valveIndex]))
                    {
                        printf("VALVE %d: PINF=%f, PSUP=%f, RATIO=%f\n",valveIndex,mean_p_inf[valveIndex],mean_p_sup[valveIndex],mean_p_inf[valveIndex]/mean_p_sup[valveIndex]);
                        printf("Crash at t=%f, ytip=%f\n",timeStepNumber*DT,y_tip);
                    }

                    // Time-average value of mean fields and mass flow rate
                    mfr[valveIndex] = compute_mfr(mean_p_inf[valveIndex],mean_p_sup[valveIndex],mean_rho_sup[valveIndex],y_tip,L_T,B0,B1,GAMMA,R);
                    
                }
            }

            /* GENERATE VALUES IN GHOST CELLS BASED ON PREVIOUS ITERATION */
            update_ghost_cells(NDOMAIN,NXtot,NYtot,x,y,xold,yold,rhoRK,uRK,vRK,pRK,HRK,ERK,TRK,B_LOC,B_TYPE,R0,NGHOST,GAMMA,R,P0,T0,M0);
            
            /* IDENTIFY SHOCK FRONTS IN EACH DIRECTION BEFORE PERFORMING AUSM-DV FLUX SPLITTING */
            for (int domainIndex = 0; domainIndex < NDOMAIN; ++domainIndex)
            {
                #pragma omp parallel for
                for (int xIndex = NGHOST; xIndex < NXtot[domainIndex]-NGHOST; ++xIndex)
                {
                    double rhoL,uL,vL,pL,rhoR,uR,vR,pR,cL,cR;

                    for (int yIndex = NGHOST; yIndex < NYtot[domainIndex]-NGHOST; ++yIndex)
                    {
                        rhoL=0.0;uL=0.0;vL=0.0;pL=0.0;rhoR=0.0;uR=0.0;vR=0.0;pR=0.0;cL=0.0;cR=0.0;

                        // MUSCL ON RIGHT FACE -> Right face flux
                        rhoL = MUSCL(rhoRK[domainIndex][xIndex-1][yIndex],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex+1][yIndex],rhoRK[domainIndex][xIndex+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[domainIndex][xIndex-1][yIndex],uRK[domainIndex][xIndex][yIndex],uRK[domainIndex][xIndex+1][yIndex],uRK[domainIndex][xIndex+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIndex][xIndex-1][yIndex],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex+1][yIndex],pRK[domainIndex][xIndex+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIndex][xIndex-1][yIndex],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex+1][yIndex],rhoRK[domainIndex][xIndex+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[domainIndex][xIndex-1][yIndex],uRK[domainIndex][xIndex][yIndex],uRK[domainIndex][xIndex+1][yIndex],uRK[domainIndex][xIndex+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIndex][xIndex-1][yIndex],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex+1][yIndex],pRK[domainIndex][xIndex+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);

                        if( ( (uL-cL)>0 && (uR-cR)<0 ) || ( (uL+cL)>0 && (uR+cR)<0 ) )
                        {
                            // A sonic point is detected in the X direction between xIndex and xIndex+1
                            #pragma omp atomic
                            sonic_x[domainIndex][xIndex][yIndex] += 1;
                            #pragma omp atomic
                            sonic_x[domainIndex][xIndex+1][yIndex] += 1;
                        }

                        // MUSCL ON LEFT FACE -> Left face flux
                        rhoL = MUSCL(rhoRK[domainIndex][xIndex-2][yIndex],rhoRK[domainIndex][xIndex-1][yIndex],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[domainIndex][xIndex-2][yIndex],uRK[domainIndex][xIndex-1][yIndex],uRK[domainIndex][xIndex][yIndex],uRK[domainIndex][xIndex+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIndex][xIndex-2][yIndex],pRK[domainIndex][xIndex-1][yIndex],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIndex][xIndex-2][yIndex],rhoRK[domainIndex][xIndex-1][yIndex],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[domainIndex][xIndex-2][yIndex],uRK[domainIndex][xIndex-1][yIndex],uRK[domainIndex][xIndex][yIndex],uRK[domainIndex][xIndex+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIndex][xIndex-2][yIndex],pRK[domainIndex][xIndex-1][yIndex],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);
                        if( ( (uL-cL)>0 && uR-cR<0 ) || ( uL+cL>0 && uR+cR<0 ) )
                        {
                            // A sonic point is detected between xIndex-1 and xIndex
                            #pragma omp atomic
                            sonic_x[domainIndex][xIndex-1][yIndex] += 1;
                            #pragma omp atomic
                            sonic_x[domainIndex][xIndex][yIndex] += 1;
                        }

                        // MUSCL ON TOP FACE -> Top face flux (Left = Down and Right = Up)
                        rhoL = MUSCL(rhoRK[domainIndex][xIndex][yIndex-1],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex][yIndex+1],rhoRK[domainIndex][xIndex][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[domainIndex][xIndex][yIndex-1],vRK[domainIndex][xIndex][yIndex],vRK[domainIndex][xIndex][yIndex+1],vRK[domainIndex][xIndex][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIndex][xIndex][yIndex-1],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex][yIndex+1],pRK[domainIndex][xIndex][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIndex][xIndex][yIndex-1],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex][yIndex+1],rhoRK[domainIndex][xIndex][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[domainIndex][xIndex][yIndex-1],vRK[domainIndex][xIndex][yIndex],vRK[domainIndex][xIndex][yIndex+1],vRK[domainIndex][xIndex][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIndex][xIndex][yIndex-1],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex][yIndex+1],pRK[domainIndex][xIndex][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);
                        if( ( (vL-cL)>0 && (vR-cR)<0 ) || ( (vL+cL)>0 && (vR+cR)<0 ) )
                        {
                            // A sonic point is detected between yIndex and yIndex+1
                            #pragma omp atomic
                            sonic_y[domainIndex][xIndex][yIndex] += 1;
                            #pragma omp atomic
                            sonic_y[domainIndex][xIndex][yIndex+1] += 1;
                        }
                        // MUSCL ON DOWN FACE -> Down face flux (Left = Down and Right = Up)
                        rhoL = MUSCL(rhoRK[domainIndex][xIndex][yIndex-2],rhoRK[domainIndex][xIndex][yIndex-1],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[domainIndex][xIndex][yIndex-2],vRK[domainIndex][xIndex][yIndex-1],vRK[domainIndex][xIndex][yIndex],vRK[domainIndex][xIndex][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIndex][xIndex][yIndex-2],pRK[domainIndex][xIndex][yIndex-1],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIndex][xIndex][yIndex-2],rhoRK[domainIndex][xIndex][yIndex-1],rhoRK[domainIndex][xIndex][yIndex],rhoRK[domainIndex][xIndex][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[domainIndex][xIndex][yIndex-2],vRK[domainIndex][xIndex][yIndex-1],vRK[domainIndex][xIndex][yIndex],vRK[domainIndex][xIndex][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIndex][xIndex][yIndex-2],pRK[domainIndex][xIndex][yIndex-1],pRK[domainIndex][xIndex][yIndex],pRK[domainIndex][xIndex][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);

                        cL = sqrt(GAMMA*pL/rhoL);
                        cR = sqrt(GAMMA*pR/rhoR);
                        if( ( (vL-cL)>0 && (vR-cR)<0 ) || ( (vL+cL)>0 && (vR+cR)<0 ) )
                        {
                            // A sonic point is detected between yIndex-1 and yIndex
                            #pragma omp atomic
                            sonic_y[domainIndex][xIndex][yIndex-1] += 1;
                            #pragma omp atomic
                            sonic_y[domainIndex][xIndex][yIndex] += 1;
                        }
                    }
                }
            }


            /* SOLVING ALL DOMAINS WITH AUSM-DV OR HANEL'S SCHEME DEPENDING ON THE POSITION OF SHOCK FRONTS */
            for (int domainIdx = 0; domainIdx < NDOMAIN; ++domainIdx)
            {
                #pragma omp parallel for
                for (int xIdx = NGHOST; xIdx < NXtot[domainIdx]-NGHOST; ++xIdx)
                {
                    double RK_k;
                    double rhoL,uL,vL,pL,rhoR,uR,vR,pR,EL,ER,HL,HR, theta;
                    double flux_l[4],flux_r[4],flux_u[4],flux_d[4],flux[4];
                    double dx,dr;
                    double y_tip, theta_tip, source[4], dx_source, dy_source;
                    int m;

                    for (int yIndex = NGHOST; yIndex < NYtot[domainIdx]-NGHOST; ++yIndex)
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
                        rhoL = MUSCL(rhoRK[domainIdx][xIdx-1][yIndex],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx+1][yIndex],rhoRK[domainIdx][xIdx+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[domainIdx][xIdx-1][yIndex],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx+1][yIndex],uRK[domainIdx][xIdx+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[domainIdx][xIdx-1][yIndex],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx+1][yIndex],vRK[domainIdx][xIdx+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIdx][xIdx-1][yIndex],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx+1][yIndex],pRK[domainIdx][xIdx+2][yIndex],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIdx][xIdx-1][yIndex],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx+1][yIndex],rhoRK[domainIdx][xIdx+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[domainIdx][xIdx-1][yIndex],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx+1][yIndex],uRK[domainIdx][xIdx+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[domainIdx][xIdx-1][yIndex],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx+1][yIndex],vRK[domainIdx][xIdx+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIdx][xIdx-1][yIndex],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx+1][yIndex],pRK[domainIdx][xIdx+2][yIndex],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // FLUX SPLITTING ON RIGHT FACE based on sonic points in Y-direction
                        if (sonic_y[domainIdx][xIdx][yIndex]+sonic_y[domainIdx][xIdx+1][yIndex]>=0)
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
                        rhoL = MUSCL(rhoRK[domainIdx][xIdx-2][yIndex],rhoRK[domainIdx][xIdx-1][yIndex],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[domainIdx][xIdx-2][yIndex],uRK[domainIdx][xIdx-1][yIndex],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[domainIdx][xIdx-2][yIndex],vRK[domainIdx][xIdx-1][yIndex],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIdx][xIdx-2][yIndex],pRK[domainIdx][xIdx-1][yIndex],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx+1][yIndex],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIdx][xIdx-2][yIndex],rhoRK[domainIdx][xIdx-1][yIndex],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[domainIdx][xIdx-2][yIndex],uRK[domainIdx][xIdx-1][yIndex],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[domainIdx][xIdx-2][yIndex],vRK[domainIdx][xIdx-1][yIndex],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIdx][xIdx-2][yIndex],pRK[domainIdx][xIdx-1][yIndex],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx+1][yIndex],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // FLUX SPLITTING ON LEFT FACE based on sonic points in Y-direction
                        if (sonic_y[domainIdx][xIdx][yIndex]+sonic_y[domainIdx][xIdx-1][yIndex]>=0)
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
                        rhoL = MUSCL(rhoRK[domainIdx][xIdx][yIndex-1],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx][yIndex+1],rhoRK[domainIdx][xIdx][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[domainIdx][xIdx][yIndex-1],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx][yIndex+1],uRK[domainIdx][xIdx][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[domainIdx][xIdx][yIndex-1],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx][yIndex+1],vRK[domainIdx][xIdx][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIdx][xIdx][yIndex-1],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx][yIndex+1],pRK[domainIdx][xIdx][yIndex+2],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIdx][xIdx][yIndex-1],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx][yIndex+1],rhoRK[domainIdx][xIdx][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[domainIdx][xIdx][yIndex-1],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx][yIndex+1],uRK[domainIdx][xIdx][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[domainIdx][xIdx][yIndex-1],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx][yIndex+1],vRK[domainIdx][xIdx][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIdx][xIdx][yIndex-1],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx][yIndex+1],pRK[domainIdx][xIdx][yIndex+2],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // PROJECT VELOCITIES IN THE FRAME OF THE UPPER FACE
                        // if (strcmp(SIM_CASE,"sup_plen_rocket")==0 || strcmp(SIM_CASE,"sup_plenum")==0)
                        // {
                        //     if (domainNumber==cone_index && xIndex==NGHOST)
                        //     {
                        //         theta = atan((y[domainNumber][xIndex][yIndex]-y[domainNumber][xIndex-1][yIndex])/(x[domainNumber][xIndex][yIndex]-x[domainNumber][xIndex-1][yIndex]));
                        //         theta_frame_velocity(uL,vL,theta);
                        //         theta_frame_velocity(uR,vR,theta);
                        //     }
                        // }

                        // FLUX SPLITTING ON TOP FACE based on sonic points in X-direction
                        if (sonic_x[domainIdx][xIdx][yIndex]+sonic_x[domainIdx][xIdx][yIndex+1]>=0)
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
                        rhoL = MUSCL(rhoRK[domainIdx][xIdx][yIndex-2],rhoRK[domainIdx][xIdx][yIndex-1],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);
                        uL = MUSCL(uRK[domainIdx][xIdx][yIndex-2],uRK[domainIdx][xIdx][yIndex-1],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);
                        vL = MUSCL(vRK[domainIdx][xIdx][yIndex-2],vRK[domainIdx][xIdx][yIndex-1],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);
                        pL = MUSCL(pRK[domainIdx][xIdx][yIndex-2],pRK[domainIdx][xIdx][yIndex-1],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx][yIndex+1],'L',MUSCL_BIAS,LIMITERNAME);

                        rhoR = MUSCL(rhoRK[domainIdx][xIdx][yIndex-2],rhoRK[domainIdx][xIdx][yIndex-1],rhoRK[domainIdx][xIdx][yIndex],rhoRK[domainIdx][xIdx][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);
                        uR = MUSCL(uRK[domainIdx][xIdx][yIndex-2],uRK[domainIdx][xIdx][yIndex-1],uRK[domainIdx][xIdx][yIndex],uRK[domainIdx][xIdx][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);
                        vR = MUSCL(vRK[domainIdx][xIdx][yIndex-2],vRK[domainIdx][xIdx][yIndex-1],vRK[domainIdx][xIdx][yIndex],vRK[domainIdx][xIdx][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);
                        pR = MUSCL(pRK[domainIdx][xIdx][yIndex-2],pRK[domainIdx][xIdx][yIndex-1],pRK[domainIdx][xIdx][yIndex],pRK[domainIdx][xIdx][yIndex+1],'R',MUSCL_BIAS,LIMITERNAME);

                        EL = pL/(GAMMA-1) + 0.5*rhoL*(pow(uL,2) + pow(vL,2));
                        ER = pR/(GAMMA-1) + 0.5*rhoR*(pow(uR,2) + pow(vR,2));
                        HL = (EL + pL)/rhoL;
                        HR = (ER + pR)/rhoR;

                        // PROJECT VELOCITIES IN THE FRAME OF THE CURRENT FACE
                        // if (strcmp(SIM_CASE,"sup_plen_rocket")==0 || strcmp(SIM_CASE,"sup_plenum")==0)
                        // {
                        //     if (domainNumber==cone_index && xIndex==NGHOST)
                        //     {
                        //         theta = atan((y[domainNumber][xIndex][yIndex-1]-y[domainNumber][xIndex-1][yIndex-1])/(x[domainNumber][xIndex][yIndex-1]-x[domainNumber][xIndex-1][yIndex-1]));
                        //         theta_frame_velocity(uL,vL,theta);
                        //         theta_frame_velocity(uR,vR,theta);
                        //     }
                        // }

                        // FLUX SPLITTING ON DOWN FACE based on sonic points in X-direction
                        if (sonic_y[domainIdx][xIdx][yIndex]+sonic_y[domainIdx][xIdx][yIndex-1]>=0)
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
                        //     if (domainNumber==cone_index && xIndex==NGHOST)
                        //     {
                        //         theta = atan((y[domainNumber][xIndex][yIndex]-y[domainNumber][xIndex-1][yIndex])/(x[domainNumber][xIndex][yIndex]-x[domainNumber][xIndex-1][yIndex]));
                        //         theta_frame_velocity(uL,vL,theta);
                        //         theta_frame_velocity(uR,vR,theta);
                        //     }
                        // }
                        // else
                        // {
                            dx = x[domainIdx][xIdx+1][yIndex] - x[domainIdx][xIdx][yIndex];
                            dr = y[domainIdx][xIdx][yIndex+1] - y[domainIdx][xIdx][yIndex];

                            flux[0] = ((flux_r[0] - flux_l[0])/dx + (flux_u[0] - flux_d[0])/dr)*DT;
                            flux[1] = ((flux_r[1] - flux_l[1])/dx + (flux_u[1] - flux_d[1])/dr)*DT;
                            flux[2] = ((flux_r[2] - flux_l[2])/dx + (flux_u[2] - flux_d[2])/dr)*DT;
                            flux[3] = ((flux_r[3] - flux_l[3])/dx + (flux_u[3] - flux_d[3])/dr)*DT;
                        // }

                        /* SOURCE TERM CORRESPONDING TO CYLINDRICAL SYSTEM OF COORDINATES */ 
                        flux[0] += rho[domainIdx][xIdx][yIndex]*v[domainIdx][xIdx][yIndex]/yc[domainIdx][xIdx][yIndex]*DT;
                        flux[1] += rho[domainIdx][xIdx][yIndex]*v[domainIdx][xIdx][yIndex]*u[domainIdx][xIdx][yIndex]/yc[domainIdx][xIdx][yIndex]*DT;
                        flux[2] += rho[domainIdx][xIdx][yIndex]*v[domainIdx][xIdx][yIndex]*v[domainIdx][xIdx][yIndex]/yc[domainIdx][xIdx][yIndex]*DT;
                        flux[3] += rho[domainIdx][xIdx][yIndex]*v[domainIdx][xIdx][yIndex]*(H[domainIdx][xIdx][yIndex])/yc[domainIdx][xIdx][yIndex]*DT;

                        /* IF NECESSARY, ADD/REMOVE FLUX AT INTERFACES */
                        if (SOLID_ON==1 && ((domainIdx==dom_low && yIndex==NYtot[dom_low]-NGHOST-1) || (domainIdx==dom_up && yIndex==NGHOST)))
                        {
                            // Locate associated valve and source term value
                            m = 0;
                            while( !(xIdx>= mfr_index_inf[m] && xIdx<=mfr_index_sup[m]) && m<N_VALVE)
                            {
                                m++;
                            }

                            // Add or remove according to domain (inverse convention because eventually substracted)
                            if (m<N_VALVE)
                            {
                                dx_source = x[domainIdx][xIdx+1][yIndex]-x[domainIdx][xIdx][yIndex];
                                dy_source = y[domainIdx][xIdx][yIndex+1]-y[domainIdx][xIdx][yIndex];
                                y_tip = U2_DOF[m][N_DOF_PER_NODE*(N_NODE-1)];
                                theta_tip = U2_DOF_K[m][N_DOF_PER_NODE*(N_NODE-1)+1];
                                compute_cell_source(domainIdx,dom_low,dom_up,r_low,r_up,mfr[m],y_tip,theta_tip,dy_source, N_V_PER_STAGE,GAMMA,R0,L_HOLE,HOLE_FACTOR,dx_source,B_HOLE,mfr_n[m],mean_rho_sup[m],mean_p_sup[m],mean_p_inf[m],rhoRK[dom_up][xIdx][yIndex],uRK[dom_up][xIdx][yIndex],vRK[dom_up][xIdx][yIndex],pRK[dom_up][xIdx][yIndex],B0,B1,L_T,source);

                                flux[0] -= source[0]*DT;
                                flux[1] -= source[1]*DT;
                                flux[2] -= source[2]*DT;
                                flux[3] -= source[3]*DT;
                                // printf("%f %f %f %f \n",source[0],source[1],source[2],source[3]);

                                if (isnan(source[0]) || isnan(source[1]) || isnan(source[2]) || isnan(source[3]))
                                {
                                    printf("DOM UP IS %d\n", dom_up);
                                    fprintf(stderr,"\nERROR:SOURCE IS NaN FOR VALVE %d at time %d and at domain %d %d %d: %f %f %f %f (ytip is %f)!\n",m,timeStepNumber,domainIdx,xIdx,yIndex,source[0],source[1],source[2],source[3],y_tip);
                                    exit(0);
                                }
                            }
                        }

                        /* INCREMENT RUNGE-KUTTA VARIABLES WITH CORRESPONDING FLUX */
                        RK_k = 1./(RK_ORDER-rungeKuttaIterationNumber);
                        rhot[domainIdx][xIdx][yIndex] = rho[domainIdx][xIdx][yIndex] - RK_k*flux[0];
                        ut[domainIdx][xIdx][yIndex] = (rho[domainIdx][xIdx][yIndex]*u[domainIdx][xIdx][yIndex] - RK_k*flux[1])/rhot[domainIdx][xIdx][yIndex];
                        vt[domainIdx][xIdx][yIndex] = (rho[domainIdx][xIdx][yIndex]*v[domainIdx][xIdx][yIndex] - RK_k*flux[2])/rhot[domainIdx][xIdx][yIndex];
                        Et[domainIdx][xIdx][yIndex] = E[domainIdx][xIdx][yIndex] - RK_k*flux[3];

                        /* DEDUCE REDUNDANT VARIABLES */
                        pt[domainIdx][xIdx][yIndex] = (GAMMA-1)*(Et[domainIdx][xIdx][yIndex] - 0.5*rhot[domainIdx][xIdx][yIndex]*(pow(ut[domainIdx][xIdx][yIndex],2)+pow(vt[domainIdx][xIdx][yIndex],2)));
                        Tt[domainIdx][xIdx][yIndex] = pt[domainIdx][xIdx][yIndex]/(R*rhot[domainIdx][xIdx][yIndex]);
                        Ht[domainIdx][xIdx][yIndex] = (Et[domainIdx][xIdx][yIndex]+pt[domainIdx][xIdx][yIndex])/rhot[domainIdx][xIdx][yIndex];
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
        // for (int xIndex = NGHOST; xIndex < NYtot[dom_low]-NGHOST; ++xIndex)
        // {
        //     p_wall += p[dom_low][NGHOST][xIndex]*(pow(y[dom_low][NGHOST][xIndex+1],2)-pow(y[dom_low][NGHOST][xIndex],2))/pow(R0,2);
        //     mfr_intake += rho[dom_low][NXtot[dom_low]-NGHOST-1][xIndex]*u[dom_low][NXtot[dom_low]-NGHOST-1][xIndex]*M_PI*(pow(y[dom_low][NXtot[dom_low]-NGHOST][xIndex+1],2)-pow(y[dom_low][NXtot[dom_low]-NGHOST][xIndex],2));
        // }
        p_wall = average_face(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],p[dom_low],NGHOST,NGHOST);
        mfr_intake = mfr_face(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],rho[dom_low],u[dom_low],NXtot[dom_low]-NGHOST-1,NGHOST);
        p_tube = domain_average(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],p[dom_low],NGHOST);
        rho_tube = domain_average(NXtot[dom_low],NYtot[dom_low],x[dom_low],y[dom_low],rho[dom_low],NGHOST);

        append_data_to_file(p_wall,DT*timeStepNumber,P_MEAN_WALL_FILENAME,EXP_EXTENSION);
        append_data_to_file(mfr_intake,DT*timeStepNumber,INTAKE_MFR_FILENAME,EXP_EXTENSION);
        append_data_to_file(p_tube,DT*timeStepNumber,P_TUBE_FILENAME,EXP_EXTENSION);
        append_data_to_file(rho_tube,DT*timeStepNumber,RHO_TUBE_FILENAME,EXP_EXTENSION);

        if (PLENUM_ON==1)
        {
            mfr_plenum = 0.0, p_plenum = 0.0, rho_plenum = 0.0, p_drag = 0.0;

            mfr_plenum = mfr_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],u[dom_up],NGHOST,NGHOST);
            p_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NGHOST);
            rho_plenum = domain_average(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],rho[dom_up],NGHOST);
            p_drag = average_face(NXtot[dom_up],NYtot[dom_up],x[dom_up],y[dom_up],p[dom_up],NXtot[dom_up]-NGHOST-1,NGHOST);

            append_data_to_file(p_plenum,DT*timeStepNumber,PLEN_P_FILENAME,EXP_EXTENSION);
            append_data_to_file(rho_plenum,DT*timeStepNumber,PLEN_RHO_FILENAME,EXP_EXTENSION);
            append_data_to_file(mfr_plenum,DT*timeStepNumber,PLENUM_MFR_FILENAME,EXP_EXTENSION);
            append_data_to_file(p_drag,DT*timeStepNumber,PLENUM_DRAG_PRESSURE_FILENAME,EXP_EXTENSION);
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
            append_data_to_file(total_mfr,DT*timeStepNumber,MFR_TOT_FILENAME,EXP_EXTENSION);
            append_data_to_file(total_pfr,DT*timeStepNumber,PFR_FILENAME,EXP_EXTENSION);
            append_multidata_to_file(N_VALVE,stage_mfr,DT*timeStepNumber,MFR_FILENAME,EXP_EXTENSION);

            // Valve displacement and pressure ratio export
            append_multidata_to_file(N_VALVE,ytip,DT*timeStepNumber,Y_TIP_FILENAME,EXP_EXTENSION);
            append_multidata_to_file(N_VALVE,pratio,DT*timeStepNumber,P_RATIO_FILENAME,EXP_EXTENSION);
        }

        /* EXPORT ALL FLUID AND SOLID DATA EVERY N STEPS */
        if(timeStepNumber%N_EXPORT==0)
        {
            // Fluid data export
            export_fluid_data(NDOMAIN,NXtot,NYtot,x,y,xc,yc,rho,u,v,p,E,T,H,DT*timeStepNumber,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT);

            // Valve data export
            if (SOLID_ON==1)
            {
                // Valve displacement
            	export_valve_data(N_VALVE,N_FEM+1,x_FEM,y_FEM,R0,DT*timeStepNumber,OUT_FOLDERNAME,EXP_EXTENSION,W_FORMAT,NGHOST);
            }
        }


        /* DISPLAY CFL */
        if(timeStepNumber%N_CFL==0)
        {
            // COMPUTE CFL (BASED ON TOTAL VELOCITY)
            CFL = get_cfl(x,y,u,v,T,R,GAMMA,NDOMAIN,NXtot,NYtot,DT);

            // COUNTING TIME BETWEEN ITERATIONS
            clock_gettime(CLOCK_REALTIME,&end_time_loop);
            SIM_TIME_LOOP = (end_time_loop.tv_sec - start_time_loop.tv_sec)+(end_time_loop.tv_nsec - start_time_loop.tv_nsec)/1E9;
            printf("Iteration %d (%.1f %%)... Max. CFL at t=%.5f sec is: %f. \n[completed in %f sec]\n",timeStepNumber,(double)(timeStepNumber)/Ntstep*100.0,DT*timeStepNumber,CFL,SIM_TIME_LOOP);
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

