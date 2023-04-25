#include "sim_case.h"
#include <stdexcept>

#include "IValve.h"
#include "parameters.h"  // List of parameters specified by user
#include "microwave.h"

SimCase::SimCase(const double simulationDuration, const double dt) :
	simulationDuration(simulationDuration),
	dt(dt)
{
	totalSimulationTimeStepCount = (int)(1 + simulationDuration / dt);  // Number of time steps
	runtimeParameters = RuntimeParameters(totalSimulationTimeStepCount);
}

void SimCase::RegisterValve(const IValve& valve)
{
	// Create a new one by copy.
	valves.push_back(valve);
	valves.back().OnRegister();
}

Domain* SimCase::AddDomain(const int id, const std::string name)
{
	//todo: Handle creation of domains, and linking them here.
	//auto it = domains.insert({ id, Domain(name) });
	domainIDS.insert({ name, id });
	return &it.first->second;
}

void SimCase::ConnectBoundaries(const int domainOneIdx, const EBoundaryLocation domainOneLocation, const int domainTwoIdx, const EBoundaryLocation domainTwoLocation)
{
	if (domainOneIdx == domainTwoIdx)
	{
		throw std::invalid_argument("Cannot connect two boundaries on the same domain");
	}

	Boundary* b1 = &domains.at(domainOneIdx).boundaries[domainOneLocation];
	Boundary* b2 = &domains.at(domainTwoIdx).boundaries[domainTwoLocation];

	// Boundaries can only be connected if they align axis-wise.
	int axis1 = (domainOneLocation < 1);
	int axis2 = (domainTwoLocation < 1);
	if (axis1 != axis2)
	{
		throw std::invalid_argument("Cannot connect two boundaries that are not aligned.");
	}

	// Can only connect if they have the same size
	if (domains.at(domainOneIdx).size[axis1] != domains.at(domainTwoIdx).size[axis1])
	{
		throw std::invalid_argument("Cannot connect two boundaries that are not the same size.");
	}

	// Can only connect if they have the same resolution
	if (domains.at(domainOneIdx).amountOfCells[axis1] != domains.at(domainTwoIdx).amountOfCells[axis1])
	{
		throw std::invalid_argument("Cannot connect two boundaries that do not have the same amount of cells.");
	}

	b1->boundaryType = EBoundaryCondition::CONNECTED;
	b2->boundaryType = EBoundaryCondition::CONNECTED;
	b1->connectedBoundary = b2;
	b2->connectedBoundary = b1;
}

void SimCase::ConnectBoundaries(const std::string domainOneName, const EBoundaryLocation domainOneLocation, const std::string domainTwoName, const EBoundaryLocation domainTwoLocation)
{
	int domainOneIdx = domainIDS.at(domainOneName);
	int domainTwoIdx = domainIDS.at(domainTwoName);
	ConnectBoundaries(domainOneIdx, domainOneLocation, domainTwoIdx, domainTwoLocation);
}

void SimCase::ApplyInitialConditions()
{
	// Setting the initial values for all the cells in the domains.
	for (auto& domainIter : domains)
	{
		Domain& domain = domainIter.second;
		switch (domain.initialisationMethod)
		{
		case EInitialisationMethod::AMBIENT_CONDITIONS:
			domain.SetToAmbientConditions(ambientTemperature, ambientStaticPressure, 0, 0, R, GAMMA);
			break;
		case EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION:
			// TODO: Move logic for determining tube length and radius to this level to allow for standing-up rockets too.
			const double tubeLength = domain.size[0];
			auto detonationConditions = SolveChapmanJougetDetonationProblem(ambientTemperature, ambientStaticPressure, ETA, S0, R, GAMMA, tubeLength, domain.size[1]);
			InitialiseDomainFromChapmanJougetDetonationSolution(&domain, detonationConditions, GAMMA);
		default:
			throw std::logic_error("The provided type of initialisation method is not (yet) implemented.");
		}
	}

	// Setting the starting deflections for the valves based on the starting fields.
	for (IValve& valve : valves)
	{
		valve.SetInitialConditions();
	}
}
