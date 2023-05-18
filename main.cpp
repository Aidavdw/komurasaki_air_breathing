/* 
* Rewrite of Florian's code, in cpp
*/

#include <ctime>
#include <iostream>

#include "sim_case.h"
#include "case_det_tube.cpp" // Temporary hard-code of specific case implementation
#include <chrono>
#include "reed_valve.h"


/* MAIN ROUTINE */
int main()
{
    /* STARTING PROGRAM, GETTING TIME... */
    std::cout << "Program starting...\n";
    std::time_t timeAtStartOfProgram = std::time(0);

    // INITIALIZE CASE
    SimCase simCase(2, 0.01);
    LoadExampleCaseWithoutReedValves(&simCase);

    /* DISPLAY INFORMATION ON TERMINAL */
    auto timeAtStartOfCalculations = std::chrono::system_clock::now();
    std::cout << "Total length of the simulation is " << simCase.simulationDuration << "s, with dt = " << simCase.dt/1000 << " ms (" << simCase.totalSimulationTimeStepCount << " steps)." << std::endl;
    //std::cout << "CFL will be displayed every " << simCase.runtimeParameters.numberOfIterationsBetweenCFLLog << "timesteps (" << simCase.runtimeParameters.numberOfIterationsBetweenCFLLog * DT*1000 <<"ms)." << std::endl;
    //std::cout << "Field properties will be exported every " << simCase.runtimeParameters.numberOfIterationsBetweenDataExport << "timesteps (" << simCase.runtimeParameters.numberOfIterationsBetweenDataExport * DT*1000 <<"ms)." << std::endl;
    
    /* INITIAL CONDITIONS ON DOMAINS */
    simCase.ApplyInitialConditionsToDomains();
    std::cout << "Initial conditions applied." << std::endl;
    
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
        for (int rungeKuttaIterationNumber = 0; rungeKuttaIterationNumber < simCase.solverSettings.rungeKuttaOrder; ++rungeKuttaIterationNumber)
        {
            // Set the ghost cells values. Note that this is deliberately done in a separate loop so that the runge-kutta operation is sure to be finished before reading into ghost cells.
            for (auto& domainIter : simCase.domains)
            {
                Domain& domain = domainIter.second;
                domain.UpdateGhostCells();
                domain.PopulateFlowDeltaBuffer(simCase.dt);
            }

            for (IValve& valve : simCase.valves)
            {
                // Update the deflections
                valve.UpdateValveState();

                // Compute mean pressure on each side and mass-flow rate at valve
                valve.CalculateFlow();
            }

            // Async await until both the buffers have been set for all the FieldQuantities in a domain. Easiest way to do this; wait until they're all finished.
            for (IValve& valve : simCase.valves)
            {
                // Add the delta due to the valve sourcing into the delta flow buffer.
                valve.AddBufferTermsToSourceCells(EFieldQuantityBuffer::FLUX);
            }
   
            for (auto& domainIter : simCase.domains)
            {
                Domain& domain = domainIter.second;
                domain.SetNextTimeStepValuesBasedOnRungeKuttaAndDeltaBuffers(rungeKuttaIterationNumber);
                domain.EmptyFlowDeltaBuffer();
            }
            
            for (IValve& valve : simCase.valves)
            {
                valve.EmptyBuffer();
            }
        } /* END OF THE RUNGE-KUTTA LOOP */
        
        // Take the 'next time step' solution of the last runge-kutta iteration as the solution for the new time step.
        for (auto& domainIter : simCase.domains)
        {
            Domain& domain = domainIter.second;
            domain.CopyFieldQuantitiesToBuffer(NEXT_TIME_STEP, CURRENT_TIME_STEP);
        }

        // todo: export flow data of domains
        // todo: export valve data
        //todo: implement displaying/calculating cfl
    } // END OF TIME LOOP


    /* FINAL INFORMATION DISPLAY AT END OF SIMULATION */
    auto timeAtFinishOfProgram = std::chrono::system_clock::now();
    // todo: implement timing and intermediate output.
}