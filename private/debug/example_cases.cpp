#include "debug/example_cases.h"
#include "domain.h"
#include "sim_case.h"

void LoadExampleCaseWithoutReedValves(SimCase* simCase) {
	// Reed valve settings
	//simCase->reed_valve_total_length = L_V + L_FIX;

	const double lengthOfTube = 0.5;
	const double heightOfTube = 0.028;
	const MeshSpacing xMeshSpacing = MeshSpacing(EMeshSpacingType::CONSTANT, lengthOfTube, 250, 0, 0);
	const MeshSpacing yMeshSpacing = MeshSpacing(EMeshSpacingType::CONSTANT, heightOfTube, 14, 0, 0);

	// Domains

	Domain* tube = simCase->AddDomain(1,
	                                  "Tube",
	                                  {0,0},
	                                  {lengthOfTube,heightOfTube},
	                                  {xMeshSpacing, yMeshSpacing },
	                                  EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION,
	                                  simCase->solverSettings.nGhost
	);
	tube->SetBoundaryType(EFace::LEFT, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);
	
	Domain* ambient = simCase->AddDomain(2,
	                                     "ambient",
	                                     {lengthOfTube,0},
	                                     {lengthOfTube,heightOfTube},
	                                     {xMeshSpacing, yMeshSpacing },
	                                     EInitialisationMethod::AMBIENT_CONDITIONS,
	                                     simCase->solverSettings.nGhost
	);
	ambient->SetBoundaryType(EFace::RIGHT, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);

	// Alternative way: simCase->ConnectBoundaries(1, EBoundaryLocation::RIGHT, 2, EBoundaryLocation::LEFT);
	simCase->ConnectBoundariesByName("Tube", EFace::RIGHT, "ambient", EFace::LEFT);
}

void LoadExampleCaseWithReedValves(SimCase* simCase)
{

	// Add An extra two domains on top.
	const double lengthOfTube = 0.5;
	const double heightOfTube = 0.028;
	const MeshSpacing xMeshSpacing = MeshSpacing(EMeshSpacingType::CONSTANT, lengthOfTube, 250, 0, 0);
	const MeshSpacing yMeshSpacing = MeshSpacing(EMeshSpacingType::CONSTANT, heightOfTube, 14, 0, 0);

	Domain* tube = simCase->AddDomain(1,
								  "Tube",
								  {0,0},
								  {lengthOfTube,heightOfTube},
								  {xMeshSpacing, yMeshSpacing },
								  EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION,
								  simCase->solverSettings.nGhost
	);
	tube->SetBoundaryType(EFace::LEFT, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);
	
	Domain* ambient = simCase->AddDomain(2,
										 "ambient",
										 {lengthOfTube,0},
										 {lengthOfTube,heightOfTube},
										 {xMeshSpacing, yMeshSpacing },
										 EInitialisationMethod::AMBIENT_CONDITIONS,
										 simCase->solverSettings.nGhost
	);
	ambient->SetBoundaryType(EFace::RIGHT, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);

	Domain* topLeft = simCase->AddDomain(3,
									  "AboveTube",
									  {0,0},
									  {lengthOfTube,heightOfTube},
									  {xMeshSpacing, yMeshSpacing },
									  EInitialisationMethod::AMBIENT_CONDITIONS,
									  simCase->solverSettings.nGhost
	);
	topLeft->SetBoundaryType(EFace::LEFT, EBoundaryCondition::SLIP);
	topLeft->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);
	topLeft->SetBoundaryType(EFace::BOTTOM, EBoundaryCondition::SLIP);
	
	Domain* topRight = simCase->AddDomain(2,
										 "AboveAmbient",
										 {lengthOfTube,0},
										 {lengthOfTube,heightOfTube},
										 {xMeshSpacing, yMeshSpacing },
										 EInitialisationMethod::AMBIENT_CONDITIONS,
										 simCase->solverSettings.nGhost
	);
	topRight->SetBoundaryType(EFace::RIGHT, EBoundaryCondition::SLIP);
	topRight->SetBoundaryType(EFace::TOP, EBoundaryCondition::SLIP);

	simCase->ConnectBoundariesByName("ambient", EFace::TOP, "aboveAmbient", EFace::BOTTOM);
	simCase->ConnectBoundariesByName("aboveTube", EFace::RIGHT, "aboveAmbient", EFace::LEFT);

	
	ReedValveGeometry reedValveGeometry;
	reedValveGeometry.beamProfile = EBeamProfile::STRAIGHT_DOUBLE_TAPERED;
	reedValveGeometry.freeLength = 25.0E-3;
	reedValveGeometry.rootThickness = 0.35E-3;
	reedValveGeometry.tipThickness = 0.15E-3;
	reedValveGeometry.rootWidth = 16.0E-3;
	reedValveGeometry.tipWidth = 12.0E-3;

	ReedValveEmpiricalParameters reedValveEmpiricalParameters;
	reedValveEmpiricalParameters.naturalFrequency = 470.0;
	reedValveEmpiricalParameters.rayleighDampingAlpha = 0;
	reedValveEmpiricalParameters.rayleighDampingBeta = 5.0E-6;
	reedValveEmpiricalParameters.dampingC1 = 5.0E-8;
	reedValveEmpiricalParameters.dampingC2 = 0.0E-8;
	reedValveEmpiricalParameters.dampingC3 = 0.0007;
	reedValveEmpiricalParameters.holeFactor = 0.9;
	

	MaterialProperties valveMaterialProperties;
	valveMaterialProperties.density = 4400.0;
	valveMaterialProperties.youngsModulus = 110.0E9;
	
	simCase->AddReedValve(tube, topLeft, TOP, lengthOfTube/3, reedValveGeometry, reedValveEmpiricalParameters, valveMaterialProperties, 6.0E-3, false);
	
	
}

