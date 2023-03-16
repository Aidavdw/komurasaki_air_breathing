#include "fem_deformation.h"

#define N_DOF_PER_NODE 2		// The total amount of degrees of freedom for a 2d beam element.

FemDeformation::FemDeformation(const int amountOfElementsToSplitBeamInto) :
	N_FEM(amountOfElementsToSplitBeamInto),
	N_NODE(amountOfElementsToSplitBeamInto + 1)	
{
	// 
	N_DOF = N_DOF_PER_NODE * N_NODE;
}