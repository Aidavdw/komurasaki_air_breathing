#include <ctime>
#include <iostream>
#include "sim_case.h"
#include "example_cases.h"
#include "main_library.h"

// Hello kinoshita! Good luck!

/* MAIN ROUTINE */
int main()
{
    /* STARTING PROGRAM, GETTING TIME... */
    std::cout << "Console Debug Program starting...\n";
    std::time_t timeAtStartOfProgram = std::time(0);

    // INITIALIZE CASE
    SimCase simCase(2, 0.01);
    LoadExampleCaseWithoutReedValves(&simCase); // Since this is the debug build, just do this very specific simple case.

    DoSimulation(simCase);
}