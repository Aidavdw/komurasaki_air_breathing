#pragma once
#include "domain.h"
#include <map>
#include "runtime_parameters.h"

struct IValve;

struct SolverSettings
{
	double AUSMSwitchBias						= 10.0;		// Bias parameter for AUSM switching function
	double MUSCLBias							= 1.0/3;	// MUSCL bias coefficient between upwind and downwind differences;
	EFluxLimiterType fluxLimiterType			= EFluxLimiterType::MIN_MOD; // The type of flux limiting that will be used in the flux splitting.
	double entropyFix							= 0.125;  // Parameter for entropy fix in AUSM-DV scheme
	int rungeKuttaOrder							= 4;		
};

struct SimCase {
	SimCase(const double simulationDuration, const double dt);
	
	RuntimeParameters runtimeParameters;
	SolverSettings solverSettings;


	std::map<int, Domain> domains; // The domains that are part of this SimCase
	std::map<std::string, int> domainIDS;
	std::vector<IValve> valves;

	double reference_mach = 0.;					// Reference mach number. I assume of the flow outside of the domain?. Formerly M0.
	double ambientTemperature = 288.15;			// Ambient temperature in Kelvin. Default is 288.15 K.
	double ambientStaticPressure = 101325;		// Ambient static pressure in pascals. Default is 101325 Pa.

	double simulationDuration;
	double dt;
	int totalSimulationTimeStepCount;			// The total amount of time steps that are in this simulation.



	// Reed valve setup
	int numberOfReedValves = 0;
	//double reed_valve_total_length = 0.; // Total length of the reed valve, with the fixed part and the flexible part combined.

	void RegisterValve(const IValve& valve);

	Domain* AddDomain(const int id, const std::string name);

	void ConnectBoundaries(const int domainOneIdx, const EFace domainOneLocation, const int domainTwoIdx, const EFace domainTwoLocation);
	void ConnectBoundaries(const std::string domainOneName, const EFace domainOneLocation, const std::string domainTwoName, const EFace domainTwoLocation);

	void ApplyInitialConditions();

};

