#include "fem_deformation.h"
#include <stdexcept>

#define N_DOF_PER_NODE 2		// The total amount of degrees of freedom for a 2d beam element.

FemDeformation::FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const EBeamProfile beamProfile, const double freeLength, const double fixedLength) :
	beamProfile(beamProfile),
	fixedNodes(amountOfFixedNodes)
{
	freeNodes = amountOfFreeSections + 1;
	amountOfNodes = freeNodes + amountOfFixedNodes;
	N_DOF = N_DOF_PER_NODE * amountOfNodes;

	CreateBeamSections();
	AssembleGlobalMassMatrix(globalMassMatrix);
	AssembleGlobalStiffnessMatrix(globalStiffnessMatrix);

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

void FemDeformation::AssembleGlobalMassMatrix(TwoDimensionalArray& matrixOut)
{
	matrixOut.Resize(N_DOF, N_DOF);

	for (const auto& beamSection : beamSections)
	{
		for (int k = 0; k < 4; ++k)
		{
			for (int j = 0; j < 4; ++j)
			{
				// superimposing all 4
				matrixOut(N_DOF_PER_NODE * beamSection.id + k, N_DOF_PER_NODE * beamSection.id + j) += beamSection.density * beamSection.massMatrix[k][j];
			}
		}
	}
}

void FemDeformation::AssembleGlobalStiffnessMatrix(TwoDimensionalArray& matrixOut)
{
	matrixOut.Resize(N_DOF, N_DOF);

	for (const auto& beamSection : beamSections)
	{
		for (int k = 0; k < 4; ++k)
		{
			for (int j = 0; j < 4; ++j)
			{
				// superimposing all 4
				matrixOut(N_DOF_PER_NODE * beamSection.id + k, N_DOF_PER_NODE * beamSection.id + j) += beamSection.youngsModulus * beamSection.stiffnessMatrix[k][j];
			}
		}
	}
}

void FemDeformation::AssembleDampingMatrix(TwoDimensionalArray& matrixOut)
{
	if (globalStiffnessMatrix.IsEmpty() || globalMassMatrix.IsEmpty())
		throw std::logic_error("Cannot assemble damping matrix, as the stiffness matrix or mass matrix has not yet been initialised.");

	for (int i = 0; i < N_DOF; ++i)
	{
		for (int j = 0; j < N_DOF; ++j)
		{
			DampingMatrix(i,j) = rayleighDampingAlpha * globalMassMatrix(i,j) + rayleighDampingBeta * globalStiffnessMatrix(i,j);
		}
	}
}




void build_damp_mat(double N_dof, double** C, double** K, double** M, double alpha, double beta)
{
	
}

void FemDeformation::AssembleNewmarkMatrix(TwoDimensionalArray& R1CholeskyOut, TwoDimensionalArray& R2Out, TwoDimensionalArray& R3Out, const double dt)
{
	R1CholeskyOut = TwoDimensionalArray(N_DOF, N_DOF);
	R2Out = TwoDimensionalArray(N_DOF, N_DOF);
	R3Out = TwoDimensionalArray(N_DOF, N_DOF);

	// We will do a cholesky decomposition + transposition on R1, so separate it for now.
	TwoDimensionalArray R1PreProcessing(N_DOF,N_DOF);
	for (int i = 0; i < N_DOF; ++i)
	{
		for (int j = 0; j < N_DOF; ++j)
		{
			R1PreProcessing(i,j) = globalMassMatrix(i,j) / dt / dt + DampingMatrix(i,j) / 2.0 / dt + globalStiffnessMatrix(i,j) / 3.0;
			R2Out(i,j) = 2.0 * globalMassMatrix(i,j) / dt / dt - globalStiffnessMatrix(i,j) / 3.0;
			R3Out(i,j)= -globalMassMatrix(i,j) / dt / dt + DampingMatrix(i,j) / 2.0 / dt - globalStiffnessMatrix(i,j) / 3.0;
		}
	}


}

