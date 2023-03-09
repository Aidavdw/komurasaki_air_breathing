#pragma once
#include "domain.h"
#include <map>


struct SimCase {
	SimCase()
	{
		
	}

	std::map<int, Domain> domains; // The domains that are part of this SimCase
	std::map<std::string, int> domainIDS;

	double reference_mach = 0.; // Reference mach number. I assume of the flow outside of the domain?. Formerly M0.

	

	// Reed valve setup
	int numberOfReedValves = 0;
	double reed_valve_total_length = 0.; // Total length of the reed valve, with the fixed part and the flexible part combined.

	Domain* AddDomain(int id, std::string name);

	void ConnectBoundaries(int domainOneIdx, EBoundaryLocation domainOneLocation, int domainTwoIdx, EBoundaryLocation domainTwoLocation);
	void ConnectBoundaries(std::string domainOneName, EBoundaryLocation domainOneLocation, std::string domainTwoName, EBoundaryLocation domainTwoLocation);

};