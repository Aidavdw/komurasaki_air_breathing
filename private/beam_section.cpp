#include "beam_section.h"
#include <cmath>
#include <stdexcept>

BeamSection::BeamSection(const double length, const std::pair<double,double> width, const std::pair<double,double> thickness, const double density, const double youngsModulus, const EBeamProfile beamProfile, const bool bCalculatePressureLoad, const int leftNodeIndex) :
	leftNodeIndex(leftNodeIndex),
	rightNodeIndex(leftNodeIndex+1),
	beamProfile(beamProfile),
	bHasPressureLoad(bCalculatePressureLoad),
	length(length),
	density(density),
	youngsModulus(youngsModulus)

{
	// idea for refactor, place left- and right variables into little structs that represent surface properties? Right now the 2 element array works too I guess.

	b[0] = width.first;
	b[1] = width.second;
	h[0] = thickness.first;
	h[1] = thickness.second;

	crossSectionalArea[0] = h[0] * b[0];
	crossSectionalArea[1] = h[1] * b[1];
	areaMomentOfInertia[0] = b[0] * std::pow(h[0], 3) / 12.0;
	areaMomentOfInertia[1] = b[1] * std::pow(h[1], 3) / 12.0;

	// Calculating the other properties; Right now, only the top surface area.
	switch (beamProfile)
	{
	case EBeamProfile::STRAIGHT_DOUBLE_TAPERED:
		topOrBottomSurfaceArea = 0.5 * (b[1] + b[0]) * length;
		break;
	default:
		throw std::logic_error("Calculating beam section properties not implemented for this beam profile.");	
	}
	
	PopulateStiffnessMatrix();
	PopulateMassMatrix();
}

void BeamSection::PopulateMassMatrix()
{
	// Populating the mass matrix
	massMatrix[0][0] = length * (10.0 * crossSectionalArea[1] + 3.0 * crossSectionalArea[0]) / 35.0;
	massMatrix[0][1] = std::pow(length, 2) * (15.0 * crossSectionalArea[1] + 7.0 * crossSectionalArea[0]) / 420.0;
	massMatrix[0][2] = 9.0 * length * (crossSectionalArea[1] + crossSectionalArea[0]) / 140.0;
	massMatrix[0][3] = -std::pow(length, 2) * (7.0 * crossSectionalArea[1] + 6.0 * crossSectionalArea[0]) / 420.0;

	massMatrix[1][0] = massMatrix[0][1];
	massMatrix[1][1] = std::pow(length, 3) * (3.0 * crossSectionalArea[1] + 5.0 * crossSectionalArea[0]) / 840.0;
	massMatrix[1][2] = -massMatrix[0][3];
	massMatrix[1][3] = -std::pow(length, 3) * (crossSectionalArea[0] + crossSectionalArea[1]) / 280.0;

	massMatrix[2][0] = massMatrix[0][2];
	massMatrix[2][1] = massMatrix[1][2];
	massMatrix[2][2] = length * (10.0 * crossSectionalArea[0] + 3.0 * crossSectionalArea[1]) / 35.0;
	massMatrix[2][3] = -std::pow(length, 2) * (15.0 * crossSectionalArea[0] + 7.0 * crossSectionalArea[1]) / 420.0;

	massMatrix[3][0] = massMatrix[0][3];
	massMatrix[3][1] = massMatrix[1][3];
	massMatrix[3][2] = massMatrix[2][3];
	massMatrix[3][3] = massMatrix[1][1];
}

void BeamSection::PopulateStiffnessMatrix()
{
	stiffnessMatrix[0][0] = 6.0 * (areaMomentOfInertia[0] + areaMomentOfInertia[1]) / std::pow(length, 3);
	stiffnessMatrix[0][1] = 2.0 / std::pow(length, 2) * (areaMomentOfInertia[0] + 2.0 * areaMomentOfInertia[1]);
	stiffnessMatrix[0][2] = -stiffnessMatrix[0][0];
	stiffnessMatrix[0][3] = 2.0 / std::pow(length, 2) * (2.0 * areaMomentOfInertia[0] + areaMomentOfInertia[1]);

	stiffnessMatrix[1][0] = stiffnessMatrix[0][1];
	stiffnessMatrix[1][1] = (areaMomentOfInertia[0] + 3.0 * areaMomentOfInertia[1]) / length;
	stiffnessMatrix[1][2] = -stiffnessMatrix[0][1];
	stiffnessMatrix[1][3] = (areaMomentOfInertia[0] + areaMomentOfInertia[1]) / length;

	stiffnessMatrix[2][0] = stiffnessMatrix[0][2];
	stiffnessMatrix[2][1] = stiffnessMatrix[1][2];
	stiffnessMatrix[2][2] = stiffnessMatrix[0][0];
	stiffnessMatrix[2][3] = -stiffnessMatrix[0][3];

	stiffnessMatrix[3][0] = stiffnessMatrix[0][3];
	stiffnessMatrix[3][1] = stiffnessMatrix[1][3];
	stiffnessMatrix[3][2] = stiffnessMatrix[2][3];
	stiffnessMatrix[3][3] = (3 * areaMomentOfInertia[0] + areaMomentOfInertia[1]) / length;
}
