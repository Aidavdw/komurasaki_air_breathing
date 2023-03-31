#pragma once
struct SimCase;

// Contains parameters related to the runtime; how the simulation is run, and reported back to the user.
struct RuntimeParameters
{
	explicit RuntimeParameters(const int totalSimulationTimeStepCount);

	int numberOfIterationsBetweenCFLLog;
	int numberOfIterationsBetweenDataExport;
};