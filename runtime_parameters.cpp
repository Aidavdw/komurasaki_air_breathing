#include "runtime_parameters.h"
#include <algorithm>
#include "parameters.h"
#include "sim_case.h"

RuntimeParameters::RuntimeParameters(const SimCase& simCase)
{
	numberOfIterationsBetweenCFLLog = std::max(1, simCase.totalSimulationTimeStepCount / N_STEP_CFL);       // Amount of iterations that are done between two displays of CFL
	numberOfIterationsBetweenDataExport = std::max(1, simCase.totalSimulationTimeStepCount / N_STEP_EXP);    // Amount of iterations that are done between two exports of data
}
