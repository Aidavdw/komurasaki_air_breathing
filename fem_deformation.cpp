#include "fem_deformation.h"
#include <stdexcept>

#define N_DOF_PER_NODE 2		// The total amount of degrees of freedom for a 2d beam element.

FemDeformation::FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const EBeamProfile beamProfile) :
	beamProfile(beamProfile),
	fixedNodes(amountOfFixedNodes)
{
	const int amountOfFreeNodes = amountOfFixedNodes + 1;
	amountOfNodes = amountOfFreeNodes + amountOfFixedNodes;

}

void FemDeformation::CreateBeamSections()
{
	if (beamSections.size() != 0)
		beamSections.clear();

	double currentNodePosX = 0;
	double width, thickness;

	for (int i = 0; i < amountOfNodes; i++)
	{
		// The 'fixed' sections are always considered at constant properties.
		if (i < fixedNodes)
		{
			width = rootWidth;
			thickness = rootThickness;
		}
		else
		{
			const double ratioCovered = (freeLength + fixedLength) / currentNodePosX;
			// Handling how the free element is created based on the selected profile.
			switch (beamProfile)
			{
			case STRAIGHTDOUBLETAPERED:
				width = rootWidth + (tipWidth - rootWidth)*ratioCovered;
				thickness = rootThickness + (tipThickness - rootThickness)*ratioCovered;
				break;
			default:
				throw std::logic_error("Determining beam properties is not implemented for this type of beam profile!");
			}
		}
		
		beamSections.emplace_back(currentNodePosX, rootWidth, rootThickness);

		// TODO: Right now, all nodes are exactly equidistant, and the same size. Make this variable?
		currentNodePosX += fixedLength / fixedNodes;

	}
}
