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
        // Set the starting runge-kutta conditions to be that of the solution of the previous time step.
        for (auto& domainIter : simCase.domains)
        {
            Domain& domain = domainIter.second;
            domain.CopyFieldQuantitiesToBuffer(EFieldQuantityBuffer::CURRENT_TIME_STEP, EFieldQuantityBuffer::RUNGE_KUTTA);
        }
        
        /* 4TH ORDER RUNGE-KUTTA PREDICTOR-CORRECTOR LOOP (iterate 4 times)*/
        for (int rungeKuttaIterationNumber = 0; rungeKuttaIterationNumber < simCase.rungeKuttaOrder; ++rungeKuttaIterationNumber)
        {
            // Set the ghost cells values. Note that this is deliberately done in a separate loop so that the runge-kutta operation is sure to be finished before reading into ghost cells.
            for (auto& domainIter : simCase.domains)
            {
                Domain& domain = domainIter.second;
                domain.UpdateGhostCells();
                domain.PopulateFlowDeltaBuffer(simCase.dt, TODO, TODO);
            }

            for (IValve& valve : simCase.valves)
            {
                // Update the deflections
                valve.Update();

                // Compute mean pressure on each side and mass-flow rate at valve
                // todo: value is currently just output. Store in IValve?
                valve.GetMassFlowRate();
                valve.PopulateValveDeltaBuffer();
            }

            // Async await until both the buffers have been set for all the FieldQuantities in a domain. Easiest way to do this; wait until they're all finished.
            for (auto& domainIter : simCase.domains)
            {
                Domain& domain = domainIter.second;
                domain.SetNextTimeStepValuesBasedOnRungeKuttaAndDeltaBuffers(rungeKuttaIterationNumber);
            }
        } /* END OF THE RUNGE-KUTTA LOOP */
        
        // Take the 'next time step' solution of the last runge-kutta iteration as the solution for the new time step.
        for (auto& domainIter : simCase.domains)
        {
            Domain& domain = domainIter.second;
            domain.CopyFieldQuantitiesToBuffer(NEXT_TIME_STEP, CURRENT_TIME_STEP);
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

