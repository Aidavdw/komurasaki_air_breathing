#include "sim_case.h"
#include <stdexcept>

#include "IValve.h"
#include "microwave.h"

SimCase::SimCase(const double simulationDuration, const double dt) :
	simulationDuration(simulationDuration),
	dt(dt)
{
	totalSimulationTimeStepCount = (int)(1 + simulationDuration / dt);  // Number of time steps
}

void SimCase::InsertValve(const IValve& valve)
{
	// Create a new one by copy.
	valves.push_back(valve);
	valves.back().OnRegister();
}

Domain* SimCase::AddDomain(const int id, const std::string name, const Position& position, const std::pair<double, double> sizeArg, const std::pair<int,int> amountOfCellsArg, const std::pair<MeshSpacing, MeshSpacing> meshSpacingArg, const EInitialisationMethod initialisationMethod, const int ghostCellDepth)
{
	// Ensure it doesn't have the same name or id
	if (domainIDS.count(name) > 0)
		throw std::logic_error("Cannot add two domains with the same name");
	if (domains.count(id) > 0)
		throw std::logic_error("Cannot add two domains with the same id");
#ifdef _DEBUG
	for (const auto& iter : domainIDS)
		if (iter.second == id)
			throw std::logic_error("Cannot add two domains with the same id");
#endif
	
	
	//todo: switch this to an emplace constructor so that it doesn't have to be copied over
	auto newDomain = Domain(name, this, position, sizeArg, amountOfCellsArg, meshSpacingArg, initialisationMethod, ghostCellDepth);
	auto it = domains.insert({id, newDomain});
	domainIDS.insert({ name, id });
	return &it.first->second;
}

void SimCase::ConnectBoundaries(const int domainOneIdx, const EFace domainOneLocation, const int domainTwoIdx, const EFace domainTwoLocation)
{
	if (domainOneIdx == domainTwoIdx)
	{
		throw std::invalid_argument("Cannot connect two boundaries on the same domain");
	}

	// Boundaries can only be connected if they align axis-wise.
	if (domainOneLocation != Opposite(domainTwoLocation))
	{
		throw std::logic_error("Right now, linking two non-opposite faces has been disabled, as this operation is not fully tested yet");
	}

	//todo: add test to see if their positions are aligned in x or y axis depending on face direction to catch more user errors.

	Boundary* b1 = &domains.at(domainOneIdx).boundaries[domainOneLocation];
	Boundary* b2 = &domains.at(domainTwoIdx).boundaries[domainTwoLocation];
	const int axis1 = (domainOneLocation == TOP || domainOneLocation == BOTTOM) ? 0 : 1;
	const int axis2 = (domainTwoLocation == TOP || domainTwoLocation == BOTTOM) ? 0 : 1;
	
	// Can only connect if they have the same size
	if (!IsCloseToZero(domains.at(domainOneIdx).size[axis1] - domains.at(domainTwoIdx).size[axis1]))
		throw std::invalid_argument("Cannot connect two boundaries that are not the same size.");

	// Can only connect if they have the same resolution
	if (domains.at(domainOneIdx).amountOfCells[axis1] != domains.at(domainTwoIdx).amountOfCells[axis1])
		throw std::invalid_argument("Cannot connect two boundaries that do not have the same amount of cells.");

	b1->boundaryType = EBoundaryCondition::CONNECTED;
	b2->boundaryType = EBoundaryCondition::CONNECTED;
	b1->connectedBoundary = b2;
	b2->connectedBoundary = b1;
}

void SimCase::ConnectBoundaries(const std::string domainOneName, const EFace domainOneLocation, const std::string domainTwoName, const EFace domainTwoLocation)
{
	int domainOneIdx = domainIDS.at(domainOneName);
	int domainTwoIdx = domainIDS.at(domainTwoName);
	ConnectBoundaries(domainOneIdx, domainOneLocation, domainTwoIdx, domainTwoLocation);
}

void SimCase::ApplyInitialConditionsToDomains()
{
	// Setting the initial values for all the cells in the domains.
	for (auto& domainIter : domains)
	{
		Domain& domain = domainIter.second;
		switch (domain.initialisationMethod)
		{
		case EInitialisationMethod::AMBIENT_CONDITIONS:
			domain.SetToAmbientConditions(ambientConditions.temperature, ambientConditions.staticPressure, ambientConditions.u, ambientConditions.v);
			break;
		case EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION:
			// TODO: Move logic for determining tube length and radius to this level to allow for standing-up rockets too.
			const double tubeLength = domain.size[0];
			auto detonationConditions = SolveChapmanJougetDetonationProblem(ambientConditions.temperature, ambientConditions.staticPressure, chapmanJougetInitialConditionParameters.ETA, chapmanJougetInitialConditionParameters.S0, chapmanJougetInitialConditionParameters.idealGasConstant, chapmanJougetInitialConditionParameters.gamma, tubeLength, domain.size[1]);
			InitialiseDomainFromChapmanJougetDetonationSolution(&domain, detonationConditions, chapmanJougetInitialConditionParameters.gamma);
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
