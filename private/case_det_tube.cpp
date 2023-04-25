// This file defines the case for a det_tube
#include "domain.h"
#include "sim_case.h"
#include "parameters.h"


void LoadCase(SimCase* simCase) {
	// Global parameters
	simCase->reference_mach;

	// Reed valve settings
	simCase->numberOfReedValves = 0;
	simCase->reed_valve_total_length = L_V + L_FIX;

	// Domains
	Domain* tube = simCase->AddDomain(1, "Tube");
	tube->XSTART = 0;
	tube->YSTART = 0;
	tube->XLENGTH = L_TUBE;
	tube->YLENGTH = R0;
	tube->GRID_RATIO_X = 1.0;
	tube->GRID_RATIO_Y = 1.0;
	tube->SetBoundaryType(EBoundaryLocation::LEFT, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EBoundaryLocation::TOP, EBoundaryCondition::SLIP);
	tube->SetBoundaryType(EBoundaryLocation::BOTTOM, EBoundaryCondition::SLIP);
	
	Domain* ambient = simCase->AddDomain(2, "Ambient");
	ambient->XSTART = L_TUBE;
	ambient->YSTART = 0;
	ambient->XLENGTH = L_TUBE;
	ambient->YLENGTH = R0;
	ambient->GRID_RATIO_X = 1.0;
	ambient->GRID_RATIO_Y = 1.0;
	ambient->SetBoundaryType(EBoundaryLocation::RIGHT, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EBoundaryLocation::TOP, EBoundaryCondition::SLIP);
	ambient->SetBoundaryType(EBoundaryLocation::BOTTOM, EBoundaryCondition::SLIP);

	// Alternative way: simCase->ConnectBoundaries(1, EBoundaryLocation::RIGHT, 2, EBoundaryLocation::LEFT);
	simCase->ConnectBoundaries("Tube", EBoundaryLocation::RIGHT, "Ambient", EBoundaryLocation::LEFT);
}

