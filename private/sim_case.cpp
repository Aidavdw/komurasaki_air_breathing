#include "sim_case.h"
#include <iostream>
#include <stdexcept>

#include "IValve.h"
#include "microwave.h"
#include "reed_valve.h"

SimCase::SimCase(const double simulationDuration, const double dt) :
	simulationDuration(simulationDuration),
	dt(dt)
{
	totalSimulationTimeStepCount = static_cast<int>(1 + simulationDuration / dt);  // Number of time steps
}

void SimCase::AddReedValve(Domain* domainThisValveFeedsInto, Domain* domainSourceFrom, const EFace boundary, const double positionAlongBoundary, const ReedValveGeometry& reedValveGeometry, const ReedValveEmpiricalParameters& reedValveEmpiricalParameters, const MaterialProperties materialProperties, const double lengthOfFixedSections, const bool bMirrored, const int amountOfFreeSections, const int amountOfFixedNodes)
{
	// SimCase::valves holds unique pointer references so it can dynamically cast them. So, create them on heap and save unique pointer.
	valves.push_back(std::make_unique<ReedValve>(domainThisValveFeedsInto, domainSourceFrom, boundary, positionAlongBoundary, reedValveGeometry , reedValveEmpiricalParameters, materialProperties, bMirrored, lengthOfFixedSections, amountOfFreeSections, amountOfFixedNodes));
	valves.back()->OnRegister();
}

void SimCase::AddReedValve(const std::string& domainIntoName, const std::string& domainOutOfName, const EFace boundary,
	const double positionAlongBoundary, const ReedValveGeometry& reedValveGeometry,
	const ReedValveEmpiricalParameters& reedValveEmpiricalParameters, const MaterialProperties materialProperties,
	const double lengthOfFixedSections, const bool bMirrored, const int amountOfFreeSections,
	const int amountOfFixedNodes)
{
	Domain* domainInto = &domains.at(domainIDS.at(domainIntoName));
	Domain* domainOutOf = &domains.at(domainIDS.at(domainOutOfName));
	AddReedValve(domainInto, domainOutOf, boundary, positionAlongBoundary, reedValveGeometry, reedValveEmpiricalParameters, materialProperties, lengthOfFixedSections, bMirrored, amountOfFreeSections, amountOfFixedNodes);
}

Domain* SimCase::AddDomain(const int id, const std::string& name, const Position& position, const std::pair<double, double> sizeArg, const std::pair<MeshSpacing, MeshSpacing>& meshSpacingArg, const EInitialisationMethod initialisationMethod)
{
	// todo: remove GhostCellDepth argument, because it's already available inside of the domain that this function is called from.
	
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

	// todo: check that the MeshSpacing size and the actual size are the same, or change SizeArg to only take in MeshSpacingArg.
	
	
	//todo: switch this to an emplace constructor so that it doesn't have to be copied over
	//auto newDomain = Domain(name, this, position, sizeArg, meshSpacingArg, initialisationMethod, ghostCellDepth);
	//auto it = domains.emplace({id, newDomain});
	auto it = domains.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple(name, this, position, sizeArg, meshSpacingArg, initialisationMethod, solverSettings.nGhost));
	domainIDS.insert({ name, id });
	return &it.first->second;
}

void SimCase::ConnectBoundariesById(const int domainOneIdx, const EFace domainOneLocation, const int domainTwoIdx, const EFace domainTwoLocation)
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
	
	// Check if they are already set- if so, throw error.
	if (domains.at(domainOneIdx).boundaries.count(domainOneLocation))
		throw std::logic_error("Cannot link boundaries, as boundary 1 has already been set.");
	if (domains.at(domainTwoIdx).boundaries.count(domainTwoLocation))
		throw std::logic_error("Cannot link boundaries, as boundary 2 has already been set.");

	//todo: add test to see if their positions are aligned in x or y axis depending on face direction to catch more user errors.
	const int axis1 = (domainOneLocation == TOP || domainOneLocation == BOTTOM) ? 0 : 1;
	const int axis2 = (domainTwoLocation == TOP || domainTwoLocation == BOTTOM) ? 0 : 1;
	
	// Can only connect if they have the same size
	if (!IsCloseToZero(domains.at(domainOneIdx).size[axis1] - domains.at(domainTwoIdx).size[axis1]))
		throw std::invalid_argument("Cannot connect two boundaries that are not the same size.");

	// Can only connect if they have the same resolution
	if (domains.at(domainOneIdx).amountOfCells[axis1] != domains.at(domainTwoIdx).amountOfCells[axis2])
		throw std::invalid_argument("Cannot connect two boundaries that do not have the same amount of cells.");

	domains.at(domainOneIdx).boundaries[domainOneLocation] = Boundary();
	auto& b1 = domains.at(domainOneIdx).boundaries.at(domainOneLocation);
	domains.at(domainTwoIdx).boundaries[domainTwoLocation] = Boundary();
	auto& b2 = domains.at(domainTwoIdx).boundaries.at(domainTwoLocation);

	b1.boundaryType = EBoundaryCondition::CONNECTED;
	b2.boundaryType = EBoundaryCondition::CONNECTED;
	b1.connectedBoundary = &b2; // Note that these pointers are invalidated upon copying!
	b2.connectedBoundary = &b1;
	b1.domain = &domains.at(domainOneIdx); // Note that these pointers are invalidated upon copying!
	b2.domain = &domains.at(domainTwoIdx);
	
}

void SimCase::ConnectBoundariesByName(const std::string& domainOneName, const EFace domainOneLocation, const std::string& domainTwoName, const EFace domainTwoLocation)
{
	int domainOneIdx = domainIDS.at(domainOneName);
	int domainTwoIdx = domainIDS.at(domainTwoName);
	ConnectBoundariesById(domainOneIdx, domainOneLocation, domainTwoIdx, domainTwoLocation);
}

void SimCase::ApplyInitialConditionsToDomainsAndValves()
{	
	// Setting the initial values for all the cells in the domains.
	for (auto& domainIter : domains)
	{
		Domain& domain = domainIter.second;
#ifdef _DEBUG
		std::cout << "Applying Initial conditions to domain '" << domain.name << std::endl;
#endif
		switch (domain.initialisationMethod)
		{
		case EInitialisationMethod::ZERO:
			// Everything is already set as zero input.
			// todo: add double-check to see if they're all 0 in debug build.
			break;
		case EInitialisationMethod::AMBIENT_CONDITIONS:
			domain.SetToAmbientConditions(ambientConditions.temperature, ambientConditions.staticPressure, ambientConditions.u, ambientConditions.v);
			break;
		case EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION:
			{
				// TODO: Move logic for determining tube length and radius to this level to allow for standing-up rockets too.
				const double tubeLength = domain.size[0];
				ChapmanJougetDetonationSolution detonationConditions = SolveChapmanJougetDetonationProblem(ambientConditions.temperature, ambientConditions.staticPressure, chapmanJougetInitialConditionParameters.energyAbsorptionCoefficient, chapmanJougetInitialConditionParameters.beamPower, chapmanJougetInitialConditionParameters.idealGasConstant, chapmanJougetInitialConditionParameters.gamma, tubeLength, domain.size[1]);
				InitialiseDomainFromChapmanJougetDetonationSolution(&domain, detonationConditions, chapmanJougetInitialConditionParameters.gamma);
				break;
			}
		case EInitialisationMethod::FROM_INPUT_RHO_P_U_V:
			{
				// If it's set from input data, it expects these fields to have been manually set already in a previous step.
				// todo: Add checker function that calls back to see if they have actually already been set externally.
				std::cout << "Using already set numpy input for initial conditions of domain '" << domain.name << "'" << std::endl;
				for (int xIdx = 0; xIdx < domain.p.currentTimeStep.nX; xIdx++)
				{
					for (int yIdx = 0; yIdx < domain.p.currentTimeStep.nY; yIdx++)
					{
						const double gamma = domain.SpecificHeatRatio();
						// shorthand refs
						const double p = domain.p.currentTimeStep.GetAt(xIdx, yIdx);
						const double u = domain.u.currentTimeStep.GetAt(xIdx, yIdx);
						const double v = domain.v.currentTimeStep.GetAt(xIdx, yIdx);
						const double rho = domain.rho.currentTimeStep.GetAt(xIdx, yIdx);
						
						const double E = p / (gamma - 1.0) + 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2));
						const double H = (E + p) / rho;

						domain.E.currentTimeStep(xIdx, yIdx) = E;
						domain.H.currentTimeStep(xIdx,yIdx) = H;
					}
				}
				break;
			}
		default:
			throw std::logic_error("The provided type of initialisation method is not (yet) implemented.");
		}
	}

	// Setting the starting deflections for the valves based on the starting fields.
	for (size_t i = 0; i < valves.size(); i++)
	{
#ifdef _DEBUG
		std::cout << "Applying Initial conditions to Valve " << i << std::endl;
#endif
		
		std::unique_ptr<IValve>& valve = valves.at(i);
		valve->SetInitialConditions();
	}

#ifdef _DEBUG
	std::cout << "Initial conditions applied." << std::endl;
#endif
	
}

void SimCase::AddRecord(const TwoDimensionalArray& src, const std::string& tag)
{
	// Verify the tag does not yet already exist
	if (twoDimensionalArrayRecords.count(tag))
		throw std::logic_error("A record with this tag has already been set");

	twoDimensionalArrayRecords[tag] = TwoDimensionalArrayRecord(&src);
}
