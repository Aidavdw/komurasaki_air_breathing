﻿#include "main_library.h"
#include <iostream>
#include "sim_case.h"
#include <chrono>
#include "reed_valve.h"

#ifdef _CREATE_DUMP_FILES
// Pls beware how cpp handles string literals with backslashes. 
#define DIR_TO_DUMP_TO std::string("C:/Users/Red Devil/Documents/git/komurasaki_air_breathing/debug_output/")
#endif

// Run a simulation, with a SimCase that has already been set up.
void DoSimulation(SimCase& simCase)
{
    std::cout << "Starting new simulation case." << std::endl;
    // Check if domains are properly initialised
    for (const auto& domain : simCase.domains)
    {
        std::cout << "Checking boundaries for domain '" << domain.second.name << "'" << std::endl;
        if (!domain.second.AreAllBoundariesSet())
            throw std::logic_error("Not all boundaries are set for domain " + domain.second.name);
    }


    /* DISPLAY INFORMATION ON TERMINAL */
    auto timeAtStartOfCalculations = std::chrono::system_clock::now();
    std::cout << "Total length of the simulation is " << simCase.simulationDuration << "s, with dt = " << simCase.dt/1000 << " ms (" << simCase.GetTotalSimulationTimeStepCount() << " steps)." << std::endl;
    //std::cout << "CFL will be displayed every " << simCase.runtimeParameters.numberOfIterationsBetweenCFLLog << "timesteps (" << simCase.runtimeParameters.numberOfIterationsBetweenCFLLog * DT*1000 <<"ms)." << std::endl;
    //std::cout << "Field properties will be exported every " << simCase.runtimeParameters.numberOfIterationsBetweenDataExport << "timesteps (" << simCase.runtimeParameters.numberOfIterationsBetweenDataExport * DT*1000 <<"ms)." << std::endl;
    
    /* INITIAL CONDITIONS ON DOMAINS */
    simCase.ApplyInitialConditionsToDomainsAndValves();
    
    std::cout << "Starting time loop." << std::endl;
    /* START TIME LOOP */
    for (int timeStepNumber = 1; timeStepNumber < simCase.GetTotalSimulationTimeStepCount(); ++timeStepNumber)
    {
        float percentageDone = static_cast<float>(timeStepNumber)/static_cast<float>(simCase.GetTotalSimulationTimeStepCount())*100.f;
        std::cout << "Advanced to time step " << timeStepNumber*simCase.dt << " [ t = " << timeStepNumber << " / " << simCase.GetTotalSimulationTimeStepCount() <<"], " << percentageDone <<"% done"  << std::endl;
        // Set the starting runge-kutta conditions to be that of the solution of the previous time step.
        for (auto& domainIter : simCase.domains)
        {
            Domain& domain = domainIter.second;
            domain.CopyFieldQuantitiesToBuffer(EFieldQuantityBuffer::CURRENT_TIME_STEP, EFieldQuantityBuffer::RUNGE_KUTTA);
        }
        
        /* 4TH ORDER RUNGE-KUTTA PREDICTOR-CORRECTOR LOOP (iterate 4 times)*/
        for (int rungeKuttaIterationNumber = 0; rungeKuttaIterationNumber < simCase.solverSettings.rungeKuttaOrder; ++rungeKuttaIterationNumber)
        {
#ifdef _DEBUG
            std::cout << "RK iter " << rungeKuttaIterationNumber << " on time step " << timeStepNumber << " [ t = " << timeStepNumber*simCase.dt << "]" << std::endl;
#endif
            // Set the ghost cells values. Note that this is deliberately done in a separate loop so that the runge-kutta operation is sure to be finished before reading into ghost cells.
            for (auto& domainIter : simCase.domains)
            {
                Domain& domain = domainIter.second;
                domain.UpdateGhostCells();
                domain.CacheEulerConservationTerms(simCase.dt);
            }


            for (std::unique_ptr<IValve>& valve : simCase.valves)
            {
#ifdef _DEBUG
                std::cout << "Updating deflections, and calculating flow deltas for valve " << valve->label << std::endl;
#endif
                // Update the deflections
                //TODO: Should this not be done in the main time loop?
                valve->UpdateValveState();

                // Compute mean pressure on each side and mass-flow rate at valve
                valve->CalculateFlow();
            }

            // Async await until both the buffers have been set for all the FieldQuantities in a domain. Easiest way to do this; wait until they're all finished.


            for (std::unique_ptr<IValve>& valve : simCase.valves)
            {
#ifdef _DEBUG
                std::cout << "Adding valve buffer terms to source cells for valve " << valve->label << std::endl;
#endif
                // Add the delta due to the valve sourcing into the delta flow buffer.
                valve->AddCachedTermsToDomainConservationEquations();
            }
            
            for (auto& domainIter : simCase.domains)
            {
                
                Domain& domain = domainIter.second;
                domain.SetNextTimeStepValuesBasedOnCachedEulerContinuities(rungeKuttaIterationNumber);
#ifdef _DEBUG
                domain.EmptyEulerConservationCache();
#endif
            }
            
            for (std::unique_ptr<IValve>& valve : simCase.valves)
            {
                valve->EmptyBuffer();
            }
        } /* END OF THE RUNGE-KUTTA LOOP */

#ifdef _DEBUG
        std::cout << "Taking the 'next time step' solution of the last runge-kutta iteration as the solution for the new time step." << std::endl;
#endif
        for (auto& domainIter : simCase.domains)
        {
            Domain& domain = domainIter.second;
            domain.CopyFieldQuantitiesToBuffer(EFieldQuantityBuffer::NEXT_TIME_STEP, EFieldQuantityBuffer::CURRENT_TIME_STEP);

#ifdef _CREATE_DUMP_FILES
            std::cout << "Dumping values into csv files" << std::endl;
            
            std::string pFilename = domain.name + "_p_time" + std::to_string(timeStepNumber) + ".csv";
            DumpPrettyWithGhostStrings(domain.p.currentTimeStep, DIR_TO_DUMP_TO + "GHOST_"+ pFilename, 4);
            DumpTwoDimensionalArrayToFileGrid(domain.p.currentTimeStep, DIR_TO_DUMP_TO + pFilename, 4);

            std::string rhoFilename = domain.name + "_rho_time" + std::to_string(timeStepNumber) + ".csv";
            DumpPrettyWithGhostStrings(domain.rho.currentTimeStep, DIR_TO_DUMP_TO + "GHOST_"+ rhoFilename, 4);
            DumpTwoDimensionalArrayToFileGrid(domain.rho.currentTimeStep, DIR_TO_DUMP_TO + rhoFilename, 4);

            std::string uFilename = domain.name + "_u_time" + std::to_string(timeStepNumber) + ".csv";
            DumpPrettyWithGhostStrings(domain.u.currentTimeStep, DIR_TO_DUMP_TO + "GHOST_"+ uFilename, 4);
            DumpTwoDimensionalArrayToFileGrid(domain.u.currentTimeStep, DIR_TO_DUMP_TO + uFilename, 4);

            std::string vFilename = domain.name + "_v_time" + std::to_string(timeStepNumber) + ".csv";
            DumpPrettyWithGhostStrings(domain.v.currentTimeStep, DIR_TO_DUMP_TO + "GHOST_"+ vFilename, 4);
            DumpTwoDimensionalArrayToFileGrid(domain.v.currentTimeStep, DIR_TO_DUMP_TO + vFilename, 4);
    #endif
        }

        // Check the Records to see what needs to be saved
        for (auto& it : simCase.twoDimensionalArrayRecords)
        {
            std::cout << "Saving record '" << it.first << "'" << std::endl;
            it.second.SaveRecord();
        }
            

        
        // todo: export valve data
        //todo: implement displaying/calculating cfl
    } // END OF TIME LOOP



    /* FINAL INFORMATION DISPLAY AT END OF SIMULATION */
    auto timeAtFinishOfProgram = std::chrono::system_clock::now();
    std::cout << "End of time loop!" << std::endl;
    std::cout << "Simulation finished \\(^-^)/" << std::endl;

#ifdef _DEBUG
    std::cout << "Saved records: " << std::endl;
    for (auto& it : simCase.twoDimensionalArrayRecords)
    {
        auto& record = it.second;
        auto& tag = it.first;
        std::cout << "\t2DArray: " << tag << "(" << record.records.size() <<" elements)" << std::endl;
    }
#endif

}
