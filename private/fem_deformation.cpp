#include "fem_deformation.h"

#include <assert.h>
#include <stdexcept>
#include <vector>

#define N_DOF_PER_NODE 2		// The total amount of degrees of freedom for a 2d beam element.

FemDeformation::FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const EBeamProfile beamProfile, const double freeLength, const double fixedLength, const double dt) :
	fixedNodes(amountOfFixedNodes),
	dt(dt),
	beamProfile_(beamProfile)
{
	freeNodes = amountOfFreeSections + 1;
	amountOfNodes = freeNodes + amountOfFixedNodes;
	N_DOF = N_DOF_PER_NODE * amountOfNodes;

	CreateBeamSections();
	nodePositionsRelativeToRoot.emplace_back(0,0);
	for (int i = 0; i < beamSections.size(); i++)
	{
		nodePositionsRelativeToRoot.emplace_back(nodePositionsRelativeToRoot[i].x + beamSections[i].length, 0, Opposite(boundary_));
	}

	AssembleGlobalMassMatrix(globalMassMatrix);
	AssembleGlobalStiffnessMatrix(globalStiffnessMatrix);
	AssembleDampingMatrix(dampingMatrix);
	
	globalStiffnessMatrixCholeskyDecomposed = CholeskyDecomposition(globalStiffnessMatrix, GetDOFVector());
	AssembleNewmarkMatrix(newmarkMatrixR1CholeskyDecomposed, newmarkMatrixR2, newmarkMatrixR3, globalStiffnessMatrixCholeskyDecomposed, dt);
}

// Empty constructor
FemDeformation::FemDeformation():
	fixedNodes(0),
	freeNodes(0),
	amountOfNodes(0),
	N_DOF(0),
	dt(0),
	beamProfile_(EBeamProfile::STRAIGHTDOUBLETAPERED),
	freeLength(0),
	fixedLength(0),
	rayleighDampingAlpha(0),
	rayleighDampingBeta(0),
	rootWidth(0),
	tipWidth(0),
	rootThickness(0),
	tipThickness(0)
{ }

void FemDeformation::UpdatePositions(const std::vector<double>& u1, const std::vector<double>& u2)
{
	// Depending on the magic number, either use u1 or u2.
	if (u2[(amountOfNodes-1)*N_DOF_PER_NODE] <= 0.012) // todo: figure out where this magic number comes from.
	{
		for (int i = 0; i < amountOfNodes; ++i)
		{
			const double deflection = u2[i*N_DOF_PER_NODE];
			// Clamp if they're going outside of the hole.
			nodePositionsRelativeToRoot[i].y = (deflection < 0) ? 0: deflection;
		}	
	}
	else
	{
		for (int i = 0; i < amountOfNodes; ++i)
		{
			// In florian's original code, no checking if < 0. Is this correct?
			nodePositionsRelativeToRoot[i].y = u1[i*N_DOF_PER_NODE];
		}
	}
}

void FemDeformation::CreateBeamSections()
{

	//todo: Propogate density and YoungsModulus here, as well as bIsFixed.
	if (!beamSections.empty())
		beamSections.clear();

	double currentNodePosX = 0;
	double length, width, thickness;

	// Pre-add the first beam section. 

	for (int i = 1; i < amountOfNodes; i++)
	{
		// The 'fixed' sections are always considered at constant properties.
		if (i < fixedNodes)
		{
			length = fixedLength / fixedNodes;
			width = rootWidth;
			thickness = rootThickness;
			BeamSection beamSection = BeamSection(length, {rootWidth, rootWidth}, {rootThickness, rootThickness}, density, youngsModulus);
			beamSections.emplace_back(beamSection);
		}
		else
		{
			const double ratioCovered = currentNodePosX / freeLength;
			// Handling how the free element is created based on the selected profile.
			switch (beamProfile_)
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
			
			BeamSection beamSection = BeamSection(length, {rootWidth, rootWidth}, {rootThickness, rootThickness}, density, youngsModulus);
			beamSections.emplace_back(beamSection);
		}

		// Advancing & resetting the values so it's easier to spot if there's a mistake.
		currentNodePosX += length;
		width = 0;
		thickness = 0;
	}
}

void FemDeformation::AssembleGlobalMassMatrix(TwoDimensionalArray& matrixOut) const
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

	#ifdef _DEBUG
	assert(matrixOut.HasDiagonalGrainsOnly(4));
	#endif
}

void FemDeformation::AssembleGlobalStiffnessMatrix(TwoDimensionalArray& matrixOut) const
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

	#ifdef _DEBUG
	assert(matrixOut.HasDiagonalGrainsOnly(4));
	#endif
}

void FemDeformation::AssembleDampingMatrix(TwoDimensionalArray& matrixOut)
{
	if (globalStiffnessMatrix.IsEmpty() || globalMassMatrix.IsEmpty())
		throw std::logic_error("Cannot assemble damping matrix, as the stiffness matrix or mass matrix has not yet been initialised.");

	for (int i = 0; i < N_DOF; ++i)
	{
		for (int j = 0; j < N_DOF; ++j)
		{
			dampingMatrix(i,j) = rayleighDampingAlpha * globalMassMatrix(i,j) + rayleighDampingBeta * globalStiffnessMatrix(i,j);
		}
	}
	
	#ifdef _DEBUG
	assert(matrixOut.HasDiagonalGrainsOnly());
	#endif
}

void FemDeformation::AssembleNewmarkMatrix(TwoDimensionalArray& R1CholeskyOut, TwoDimensionalArray& R2Out, TwoDimensionalArray& R3Out, TwoDimensionalArray& KCholeskyOut, const double dt) const
{
	if (globalStiffnessMatrix.IsEmpty() || globalMassMatrix.IsEmpty() || dampingMatrix.IsEmpty())
		throw std::logic_error("Cannot assemble Newmark matrix, as the stiffness matrix or mass matrix has not yet been initialised.");

	R2Out.Resize(N_DOF, N_DOF);
	R3Out.Resize(N_DOF, N_DOF);

	// We will do a cholesky decomposition + transposition on R1, so separate it for now.
	TwoDimensionalArray R1PreProcessing(N_DOF,N_DOF);
	for (int i = 0; i < N_DOF; ++i)
	{
		for (int j = 0; j < N_DOF; ++j)
		{
			R1PreProcessing(i,j) = globalMassMatrix.GetAt(i,j) / dt / dt + dampingMatrix.GetAt(i,j) / 2.0 / dt + globalStiffnessMatrix.GetAt(i,j) / 3.0;
			R2Out(i,j) = 2.0 * globalMassMatrix.GetAt(i,j) / dt / dt - globalStiffnessMatrix.GetAt(i,j) / 3.0;
			R3Out(i,j)= -globalMassMatrix.GetAt(i,j) / dt / dt + dampingMatrix.GetAt(i,j) / 2.0 / dt - globalStiffnessMatrix.GetAt(i,j) / 3.0;
		}
	}

	const auto DOFVector = GetDOFVector();
	R1CholeskyOut = CholeskyDecomposition(R1PreProcessing, DOFVector).Transpose();
	KCholeskyOut = CholeskyDecomposition(globalStiffnessMatrix, DOFVector).Transpose();

	#ifdef _DEBUG
	assert(R1CholeskyOut.IsLowerTriangular());
	assert(KCholeskyOut.IsLowerTriangular());
	#endif
}

std::vector<double> FemDeformation::GetDOFVector() const
{
	int dofIndex = 0;
	std::vector<double> DOFVector(freeNodes, 0);
	for (int i = fixedNodes; i < amountOfNodes; ++i)
	{
		for (int j = 0; j < N_DOF_PER_NODE; ++j)
		{
			DOFVector[dofIndex] = N_DOF_PER_NODE * i + j;
			dofIndex++;
		}
	}

	return DOFVector;
}

/* Given a symmetric positive definite matrix M (size NxN) (should be checked by user), this function computes the Cholesky decomposition of this matrix and fills a matrix L (that should be allocated and initialized by the user) so that A = L*LT (LT is the transpose matrix of L). The output L matrix is a superior triangular matrix. */
TwoDimensionalArray FemDeformation::CholeskyDecomposition(const TwoDimensionalArray& matrix, const std::vector<double>& DOFVector)
{
	#ifdef _DEBUG
	// Technically, the matrix m ust be positive-definite, which is a subset of symmetric matrices where all the local determinants are positive. However, testing that is relatively costly for a big matrix and its a big long algorithm, so it's left out- this is just to catch the obvious cases anyway!
	if (!matrix.IsDiagonallySymmetric())
		throw std::logic_error("An matrix that is not diagonally symmetric cannot be cholesky decomposed.");
	#endif
	
	// This is a variation on the Cholesky-Crout algorithm?
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
			throw std::logic_error("While Cholesky decomposing, got a non-positive definite matrix!");
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

	#ifdef _DEBUG
	assert(out.IsUpperTriangular());
	#endif

	// Todo: Change to return by reference, not value.
	return out;
}
void FemDeformation::SolveCholeskySystem(std::vector<double> &deflectionVectorOut, const std::vector<double>& load) const
{
	#ifdef _DEBUG
	// Check if the vector sizes are compatible
	if (static_cast<int>(load.size()) != globalStiffnessMatrix.nX)
	{
		throw std::invalid_argument("Cannot solve cholesky system, as the load vector is not the right size.");
	}

	// Check if the matrices are properly triangular.
	assert(globalStiffnessMatrixCholeskyDecomposed.IsLowerTriangular());
	#endif
	
	
	deflectionVectorOut.resize(load.size());
	std::fill(deflectionVectorOut.begin(), deflectionVectorOut.end(), 0);
	std::vector<double> DOFVector = GetDOFVector();


	// First step : Forward substitution starting from first index
	for (int i = 0; i < freeNodes * N_DOF_PER_NODE; ++i)
	{
		double sum = 0.0;
		for (int j = 0; j < i; ++j)
		{
			//TODO: check if this does not weirdly cut the matrix.
			sum += globalStiffnessMatrixCholeskyDecomposed.GetAt(i,j)*deflectionVectorOut[DOFVector.at(j)];
		}
		deflectionVectorOut[DOFVector.at(i)] = (load[DOFVector.at(i)] - sum)/globalStiffnessMatrixCholeskyDecomposed.GetAt(i,i);
	}

	
	// Second step : Backward subsitution starting from last index
	for (int i = freeNodes * N_DOF_PER_NODE - 1; i > -1; --i)
	{
		double sum = 0.0;
		for (int j = i + 1; j < freeNodes * N_DOF_PER_NODE; ++j)
		{
			sum += globalStiffnessMatrixCholeskyDecomposed.GetAt(j,i)*deflectionVectorOut[DOFVector.at(j)];
		}
		deflectionVectorOut[DOFVector.at(i)] = (deflectionVectorOut[DOFVector.at(i)] - sum)globalStiffnessMatrixCholeskyDecomposed.GetAt(i,i);
	}

	#ifdef _DEBUG
	for (const auto num : deflectionVectorOut)
		if(isnan(num))
			throw std::logic_error("Cholesky solution contains an NAN");
	#endif
}
BeamSection* FemDeformation::BeamSectionsConnectedToNode(const int nodeIndex, const bool bRight)
{
	if (nodeIndex == 0 && bRight == false)
		return nullptr;
	else if (nodeIndex == (amountOfNodes - 1) && bRight == true )
		return nullptr;
	else if (nodeIndex >= amountOfNodes)
	{
		#ifdef _DEBUG
		throw std::invalid_argument("Cannot get a beam section for a node index that is longer than the amount of nodes!");
		#endif
		return nullptr;
	}

	return &(beamSections[nodeIndex - fixedNodes]);
}