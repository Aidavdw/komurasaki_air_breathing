#pragma once
#include "domain.h"
#include <map>

#include "beam_section.h"
#include "IValve.h"
#include "reed_valve_geometry.h"
#include "pythoninterface/record.h"

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
	double beamPower          = 2000.0E3;   // Total microwave beam power (W)
	double energyAbsorptionCoefficient         = 1.0;        // Energy absorption coefficient

	// Note that these are fixed, not calculated.
	double gamma	   = 1.4;
	double idealGasConstant = 287.0;
};

// Contains parameters related to the runtime; how the simulation is run, and reported back to the user.
struct RuntimeParameters
{
	int numberOfTimeStepsBetweenCflLog = 10; 
	int numberOfTimeStepsBetweenDataExport = 10; // todo: actually respect this, and not dump data EVERY time step.
};

// Master container for configuration settings for the simulation, domains, and valves.
class SimCase {
public:
	SimCase(const double simulationDuration, const double dt);
	
	RuntimeParameters runtimeParameters;
	SolverSettings solverSettings;
	AmbientConditions ambientConditions;
	ChapmanJougetInitialConditionParameters chapmanJougetInitialConditionParameters;

	std::map<int, Domain> domains;			// The domains that are part of this SimCase. Add more using AddDomain()
	std::map<std::string, int> domainIDS;	// names of domains to their IDs.
	std::vector<std::unique_ptr<IValve>> valves;				// The valves that are part of this simulation. Add more using InsertValve-like functions, defind at the bottom. Holds unique pointers to the base class, so that they can be memory-tracked.
	std::map<std::string, TwoDimensionalArrayRecord> twoDimensionalArrayRecords;	// Contains all the records, indicating which values are tracked and saved so that they can be analysed after running the simulation.

	double simulationDuration;				// Total amount of time to run the simulation for
	double dt;								// The amount of time between two time steps.
	int totalSimulationTimeStepCount;		// The total amount of time steps that are in this simulation, = simulationDuration/dt



	Domain* AddDomain(const int id, const std::string& name, const Position& position, const std::pair<MeshSpacing, MeshSpacing>& meshSpacingArg, const EInitialisationMethod initialisationMethod); // Adds a new domain to simulation.

	void ConnectBoundariesById(const int domainOneIdx, const EFace domainOneLocation, const int domainTwoIdx, const EFace domainTwoLocation); // Sets two boundaries on two different domains to be connected, meaning that flow going out of one goes into the other (and vice versa).
	void ConnectBoundariesByName(const std::string& domainOneName, const EFace domainOneLocation, const std::string& domainTwoName, const EFace domainTwoLocation); // Sets two boundaries on two different domains to be connected, meaning that flow going out of one goes into the other (and vice versa).
	
	void ApplyInitialConditionsToDomainsAndValves(); // Sets initial conditions for all domains based on their initialisation method, then puts the reed valves in a neutral position inside of them.

	void AddRecord(const TwoDimensionalArray& src, const std::string& tag);	// Add a record to a specific variable, meaning that it will be saved for every time step so that it can be analysed after the simulation finishes.
	Domain& GetDomainByID(const int id) {return domains.at(id);}
	Domain& GetDomainByName(const std::string& name) {return GetDomainByID(domainIDS.at(name));}

	/**** Valve inserters. If you derive a new Valve class that derives from IValve, add a new function here to register it. ****/
	void AddReedValve(Domain* domainThisValveFeedsInto, Domain* domainSourceFrom, const EFace boundary, const double positionAlongBoundary, const
	                  ReedValveGeometry& reedValveGeometry, const ReedValveEmpiricalParameters& reedValveEmpiricalParameters, const
	                  MaterialProperties materialProperties, const double
	                  lengthOfFixedSections, const bool bMirrored, const int amountOfFreeSections = 30, const int amountOfFixedNodes = 3);

	void AddReedValve(const std::string& domainIntoName, const std::string& domainOutOfName, const EFace boundary, const double positionAlongBoundary, const ReedValveGeometry& reedValveGeometry, const ReedValveEmpiricalParameters& reedValveEmpiricalParameters, const MaterialProperties materialProperties, const double lengthOfFixedSections, const bool bMirrored, const int amountOfFreeSections = 30, const int amountOfFixedNodes = 3);
	
	
	 
};

