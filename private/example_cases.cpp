#include "example_cases.h"
#include "domain.h"
#include "sim_case.h"

// Thisdefines the case for a det_tube
void LoadExampleCaseWithoutReedValves(SimCase* simCase) {
	// Reed valve settings
	//simCase->reed_valve_total_length = L_V + L_FIX;

	const double lengthOfTube = 0.5;
	const double heightOfTube = 0.028;
	MeshSpacing constantMeshSpacing = MeshSpacing(EMeshSpacingType::CONSTANT, lengthOfTube, 6, 0, 0);

	// Domains

	Domain* tube = simCase->AddDomain(1,
	                                  "Tube",
	                                  {lengthOfTube/2,heightOfTube/2},
	                                  {lengthOfTube,heightOfTube},
	                                  {constantMeshSpacing, constantMeshSpacing },
	                                  EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION,
	                                  simCase->solverSettings.nGhost
	);
	tube->SetBoundaryType(EFace::LEFT, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);
	
	Domain* ambient = simCase->AddDomain(2,
	                                     "ambient",
	                                     {lengthOfTube/2,heightOfTube/2},
	                                     {lengthOfTube,heightOfTube},
	                                     {constantMeshSpacing, constantMeshSpacing },
	                                     EInitialisationMethod::AMBIENT_CONDITIONS,
	                                     simCase->solverSettings.nGhost
	);
	ambient->SetBoundaryType(EFace::RIGHT, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);

	// Alternative way: simCase->ConnectBoundaries(1, EBoundaryLocation::RIGHT, 2, EBoundaryLocation::LEFT);
	simCase->ConnectBoundaries("Tube", EFace::RIGHT, "ambient", EFace::LEFT);
}

