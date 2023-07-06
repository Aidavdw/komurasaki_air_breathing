#include <ctime>
#include <iostream>
#include "sim_case.h"
#include "debug/example_cases.h"
#include "main_library.h"

/* Note that this is **NOT** what is called by the python interface. This is what is called for the debug build. To see what gets called with a case set-up with the python interface, see DoSimulation() */
int main()
{
    /* STARTING PROGRAM, GETTING TIME... */
    std::cout << "Console Debug Program starting...\n";
    std::time_t timeAtStartOfProgram = std::time(0);

    // INITIALIZE CASE
    SimCase simCase(2, 0.01);
    LoadExampleCaseWithoutReedValves(&simCase); // Since this is the debug build, just do this very specific simple case.
    //LoadExampleCaseWithReedValves(&simCase);
    
    DoSimulation(simCase);

    std::cout << "Done!!!!" << std::endl;
}