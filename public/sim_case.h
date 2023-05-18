#pragma once
#include "domain.h"
#include <map>

// forward declarations
class IValve;

// Changes the behaviour of how the solver works.
struct SolverSettings
{
	double AUSMSwitchBias						= 10.0;		// Bias parameter for AUSM switching function
	double MUSCLBias							= 1.0/3;	// MUSCL bias coefficient between upwind and downwind differences;
	EFluxLimiterType fluxLimiterType			= EFluxLimiterType::MIN_MOD; // The type of flux limiting that will be used in the flux splitting.
	double entropyFix							= 0.125;  // Parameter for entropy fix in AUSM-DV scheme
	int rungeKuttaOrder							= 4;
	int nGhost							= 2;
};

// Represents the conditions far outside of the domain. Is used to initialise the system.
struct AmbientConditions
{
	double mach = 0.;					// Reference mach number. I assume of the flow outside of the domain?. Formerly M0.
	double temperature = 288.15;			// Ambient temperature in Kelvin. Default is 288.15 K.
	double staticPressure = 101325;		// Ambient static pressure in pascals. Default is 101325 Pa.
	double v = 0;
	double u = 0;
};

// The vales that are used if an initial condition is calculated using Chapman Jouget theory.
struct ChapmanJougetInitialConditionParameters
{
	double S0          = 2000.0E3;   // Total microwave beam power (W)
	double ETA         = 1.0;        // Energy absorption coefficient

	// Note that these are fixed, not calculated.
	double gamma	   = 1.4;
	double idealGasConstant = 287.0;
};

// Contains parameters related to the runtime; how the simulation is run, and reported back to the user.
struct RuntimeParameters
{
	int numberOfTimeStepsBetweenCflLog = 10;
	int numberOfTimeStepsBetweenDataExport = 10;
};

struct SimCase {
	SimCase(const double simulationDuration, const double dt);
	
	RuntimeParameters runtimeParameters;
	SolverSettings solverSettings;
	AmbientConditions ambientConditions;
	ChapmanJougetInitialConditionParameters chapmanJougetInitialConditionParameters;

	std::map<int, Domain> domains; // The domains that are part of this SimCase
	std::map<std::string, int> domainIDS;
	std::vector<IValve> valves;

	double simulationDuration;
	double dt;
	int totalSimulationTimeStepCount;			// The total amount of time steps that are in this simulation.

	void InsertValve(const IValve& valve);

	Domain* AddDomain(const int id, const std::string name, const Position& position, const std::pair<double, double> sizeArg, const std::
	                  pair<int, int> amountOfCellsArg, const std::pair<MeshSpacing, MeshSpacing> meshSpacingArg, const EInitialisationMethod
	                  initialisationMethod, const int ghostCellDepth);

	void ConnectBoundaries(const int domainOneIdx, const EFace domainOneLocation, const int domainTwoIdx, const EFace domainTwoLocation);
	void ConnectBoundaries(const std::string domainOneName, const EFace domainOneLocation, const std::string domainTwoName, const EFace domainTwoLocation);

	//todo: add proxy inserter for (reed) valves so that it doesn't need to be constructed so awkwardly.

	void ApplyInitialConditionsToDomains();

};

