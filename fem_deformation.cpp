#include "fem_deformation.h"
#include <stdexcept>
#include <vector>

#define N_DOF_PER_NODE 2		// The total amount of degrees of freedom for a 2d beam element.

FemDeformation::FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const EBeamProfile beamProfile, const double freeLength, const double fixedLength, const double dt) :
	beamProfile(beamProfile),
	fixedNodes(amountOfFixedNodes),
	dt(dt)
{
	freeNodes = amountOfFreeSections + 1;
	amountOfNodes = freeNodes + amountOfFixedNodes;
	N_DOF = N_DOF_PER_NODE * amountOfNodes;

	CreateBeamSections();
	nodePositionsRelativeToRoot.push_back({ 0,0 });
	for (int i = 0; i < beamSections.size(); i++)
	{
		nodePositionsRelativeToRoot.push_back({ nodePositionsRelativeToRoot[i].x + beamSections[i].length, 0 });
	}

	AssembleGlobalMassMatrix(globalMassMatrix);
	AssembleGlobalStiffnessMatrix(globalStiffnessMatrix);
	AssembleNewmarkMatrix(newmarkMatrixR1CholeskyDecomposed, newmarkMatrixR2, newmarkMatrixR3, globalStiffnessMatrixCholeskyDecomposed, dt);

}

void FemDeformation::CreateBeamSections()
{
	if (beamSections.size() != 0)
		beamSections.clear();

	double currentNodePosX = 0;
	double length, width, thickness;

	for (int i = 0; i < amountOfNodes; i++)
	{
		// The 'fixed' sections are always considered at constant properties.
		if (i < fixedNodes)
		{
			length = fixedLength / fixedNodes;
			width = rootWidth;
			thickness = rootThickness;
		}
		else
		{
			const double ratioCovered = currentNodePosX / freeLength;
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
			// TODO: Right now, all nodes are exactly equidistant, and the same size. Make this variable?
			length = freeLength / freeNodes;
		}
		
		beamSections.emplace_back(length, rootWidth, rootThickness);

		currentNodePosX += length;

	}
}

void FemDeformation::AssembleGlobalMassMatrix(TwoDimensionalArray& matrixOut)
{
	matrixOut.Resize(N_DOF, N_DOF);

	for (int beamIdx = 0; beamIdx < beamSections.size(); beamIdx++)
	{
		const auto& beamSection = beamSections[beamIdx];
		for (int k = 0; k < 4; ++k)
		{
			for (int j = 0; j < 4; ++j)
			{
				// superimposing all 4
				matrixOut(N_DOF_PER_NODE * beamIdx + k, N_DOF_PER_NODE * beamIdx + j) += beamSection.density * beamSection.massMatrix[k][j];
			}
		}
	}
}

void FemDeformation::AssembleGlobalStiffnessMatrix(TwoDimensionalArray& matrixOut)
{
	matrixOut.Resize(N_DOF, N_DOF);

	for (int beamIdx = 0; beamIdx < beamSections.size(); beamIdx++)
	{
		const auto& beamSection = beamSections[beamIdx];
		for (int k = 0; k < 4; ++k)
		{
			for (int j = 0; j < 4; ++j)
			{
				// superimposing all 4
				matrixOut(N_DOF_PER_NODE * beamIdx + k, N_DOF_PER_NODE * beamIdx + j) += beamSection.youngsModulus * beamSection.stiffnessMatrix[k][j];
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

void FemDeformation::AssembleNewmarkMatrix(TwoDimensionalArray& R1CholeskyOut, TwoDimensionalArray& R2Out, TwoDimensionalArray& R3Out, TwoDimensionalArray& KCholeskyOut, const double dt)
{
	if (globalStiffnessMatrix.IsEmpty() || globalMassMatrix.IsEmpty() || DampingMatrix.IsEmpty())
		throw std::logic_error("Cannot assemble newmark matrix, as the stiffness matrix or mass matrix has not yet been initialised.");

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

	auto DOFVector = GetDOFVector();
	R1CholeskyOut = CholeskyDecomposition(R1PreProcessing, DOFVector).Transpose();
	KCholeskyOut = CholeskyDecomposition(globalStiffnessMatrix, DOFVector).Transpose();
}

std::vector<double> FemDeformation::GetDOFVector()
{
	int dofIndex = 0;
	std::vector<double> dof_vector(freeNodes, 0);
	for (int i = fixedNodes; i < amountOfNodes; ++i)
	{
		for (int j = 0; j < N_DOF_PER_NODE; ++j)
		{
			dof_vector[dofIndex] = N_DOF_PER_NODE * i + j;
			dofIndex++;
		}
	}
}

/* Given a symmetric positive definite matrix M (size NxN) (should be checked by user), this function computes the Cholesky decomposition of this matrix and fills a matrix L (that should be allocated and initialized by the user) so that A = L*LT (LT is the transpose matrix of L). The outputed L matrix is a superior triangular matrix. */
TwoDimensionalArray CholeskyDecomposition(const TwoDimensionalArray& matrix, std::vector<double>& DOFVector)
{
	//const int activeDegreesOfFreedom = freeNodes * N_DOF_PER_NODE; <- the original thing, but should work like this too.
	const int activeDegreesOfFreedom = matrix.nX;
	TwoDimensionalArray out(matrix.nX, matrix.nX);
	for (int i = 0; i < activeDegreesOfFreedom; ++i)
	{
		out(i,i) = matrix.GetAt(DOFVector[i],DOFVector[i]);
		for (int k = 0; k < i; ++k)
		{
			out(i,i) -= out(k,i) * out(k,i);
		}
		if (out(i,i) <= 0)
		{
			throw std::logic_error("While cholesky decomposing, got a non-positive definite matrix!");
		}
		out(i,i) = sqrt(out(i, i));

		for (int j = i + 1; j < activeDegreesOfFreedom; ++j)
		{
			out(i, j) = matrix.GetAt(DOFVector[i], DOFVector[j]);
			for (int k = 0; k < i; ++k)
			{
				out(i,j) -= out(k,i) * out(k,j);
			}
			out(i,j) /= out(i,i);
		}
	}

	return out;
}