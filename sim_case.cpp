#include "sim_case.h"
#include <stdexcept>

#include "parameters.h"  // List of parameters specified by user

SimCase::SimCase()
{
	totalSimulationTimeStepCount = (int)(1 + TSIM / DT);  // Number of time steps
	runtimeParameters = RuntimeParameters(this);

}

Domain* SimCase::AddDomain(const int id, const std::string name)
{
	auto it = domains.insert({ id, Domain(name) });
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

	b1->boundaryType = EBoundaryType::CONNECTED;
	b2->boundaryType = EBoundaryType::CONNECTED;
	b1->connectedBoundary = b2;
	b2->connectedBoundary = b1;
}

void SimCase::ConnectBoundaries(const std::string domainOneName, const EBoundaryLocation domainOneLocation, const std::string domainTwoName, const EBoundaryLocation domainTwoLocation)
{
	int domainOneIdx = domainIDS.at(domainOneName);
	int domainTwoIdx = domainIDS.at(domainTwoName);
	ConnectBoundaries(domainOneIdx, domainOneLocation, domainTwoIdx, domainTwoLocation);
}
