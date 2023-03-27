#pragma once
#include "domain.h"
#include <map>
#include "runtime_parameters.h"


struct SimCase {
	SimCase();
	
	RuntimeParameters runtimeParameters;

	std::map<int, Domain> domains; // The domains that are part of this SimCase
	std::map<std::string, int> domainIDS;

	double reference_mach = 0.;					// Reference mach number. I assume of the flow outside of the domain?. Formerly M0.
	double ambientTemperature = 288.15;			// Ambient temperature in Kelvin. Default is 288.15 K.
	double ambientStaticPressure = 101325;		// Ambient static pressure in pascals. Default is 101325 Pa.


	int totalSimulationTimeStepCount;			// The total amount of time steps that are in this simulation.

	// Reed valve setup
	int numberOfReedValves = 0;
	double reed_valve_total_length = 0.; // Total length of the reed valve, with the fixed part and the flexible part combined.

	void RegisterValves(const std::vector<Valve>& valves);

	Domain* AddDomain(const int id, const std::string name);

	void ConnectBoundaries(const int domainOneIdx, const EBoundaryLocation domainOneLocation, const int domainTwoIdx, const EBoundaryLocation domainTwoLocation);
	void ConnectBoundaries(const std::string domainOneName, const EBoundaryLocation domainOneLocation, const std::string domainTwoName, const EBoundaryLocation domainTwoLocation);

	void ApplyInitialConditions();

};