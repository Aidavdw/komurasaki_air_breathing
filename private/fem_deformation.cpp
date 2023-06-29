#include "fem_deformation.h"

#include <cassert>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>

#include "AuxFunctions.h"

FemDeformation::FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const double fixedLength, const ReedValveGeometry& reedValveGeometry , const ReedValveEmpiricalParameters& reedValveEmpiricalParameters, const MaterialProperties& materialProperties, const double dt, const EDirection upDirection) :
	fixedNodes(amountOfFixedNodes),
	geometry(reedValveGeometry),
	empiricalParameters(reedValveEmpiricalParameters),
	material(materialProperties),
	dt(dt),
	rootWidth(reedValveGeometry.rootWidth),
	tipWidth(reedValveGeometry.tipWidth),
	rootThickness(reedValveGeometry.rootThickness),
	tipThickness(reedValveGeometry.tipThickness),
	freeLength(reedValveGeometry.freeLength),
	fixedLength(fixedLength)
{
	freeNodes = amountOfFreeSections + 1;
	amountOfNodes = freeNodes + amountOfFixedNodes;
	N_DOF = N_DOF_PER_NODE * amountOfNodes;
	n_active = freeNodes*N_DOF_PER_NODE;

	CreateBeamSections();
	nodePositionsRelativeToRoot.emplace_back(0,0);
	for (int i = 0; i < beamSections.size(); i++)
	{
		//todo: check if upDirection should not be mirrored, due to the new coordinate way of handling
		nodePositionsRelativeToRoot.emplace_back(nodePositionsRelativeToRoot[i].x + beamSections[i].length, 0, upDirection);
	}

	AssembleGlobalMassMatrix(globalMassMatrix);
	AssembleGlobalStiffnessMatrix(globalStiffnessMatrix);
	AssembleDampingMatrix(dampingMatrix);
	
	globalStiffnessMatrixCholeskyDecomposed = CholeskyDecomposition(globalStiffnessMatrix, GetDOFVector());
	AssembleNewmarkMatrix(newmarkMatrixR1CholeskyDecomposed, newmarkMatrixR2, newmarkMatrixR3, globalStiffnessMatrixCholeskyDecomposed, dt);

	positionsInPreviousTimeStep = std::vector<Position>(amountOfNodes);
	nodePositionsRelativeToRoot = std::vector<Position>(amountOfNodes);
}

// Empty constructor
FemDeformation::FemDeformation() :
	fixedNodes(0),
	freeNodes(0),
	amountOfNodes(0),
	N_DOF(0),
	n_active(0),
	dt(0),
	rootWidth(0),
	tipWidth(0),
	rootThickness(0),
	tipThickness(0),
	freeLength(0),
	fixedLength(0)
{ }

void FemDeformation::UpdatePositions(const std::vector<double>& newDeflection)
{
	#ifdef _DEBUG
	assert(nodePositionsRelativeToRoot.size() == positionsInPreviousTimeStep.size());
	assert(nodePositionsRelativeToRoot.size() == 2*newDeflection.size());
	#endif
	// First set the values of the previous time step as the one from the current time step.
	std::copy(nodePositionsRelativeToRoot.begin(), nodePositionsRelativeToRoot.end(), std::back_inserter(positionsInPreviousTimeStep));

	#ifdef _DEBUG
	for (int i = 0; i<nodePositionsRelativeToRoot.size(); i++)
	{
		assert(IsCloseToZero(nodePositionsRelativeToRoot[i].x - positionsInPreviousTimeStep[i].x));
		assert(IsCloseToZero(nodePositionsRelativeToRoot[i].y - positionsInPreviousTimeStep[i].y));
	}
	#endif
	
	// in the deflection vector, it is x1,y1,x2,x3,[...], so skip an index every time
	for (int i = 0; i < nodePositionsRelativeToRoot.size(); i++)
	{
		double nodeYDeflection = newDeflection[i*2 + 1];
		nodeYDeflection = (nodeYDeflection < 0) ? 0: nodeYDeflection;			//clamp to be minimum 0, don't allow going outside of the hole backwards.
		nodeYDeflection = (nodeYDeflection > 0.012) ? 0.012: nodeYDeflection;	// clamp to max 0.012, magic value. Not sure where it comes from.
		nodePositionsRelativeToRoot[i].y = newDeflection[i*2 + 1];
		
		// Not sure if I should set the x positions too, as florian's code does NOT do this.
		//nodePositionsRelativeToRoot[i].x = newDeflection[i*2];
	}

}
void FemDeformation::CalculateNewDeflections(std::vector<double> &u2Out, const std::vector<double> &load) const
{
	std::vector<double> deflectionNow, deflectionPrevious;
	GetDeflectionVectorFromPositions(deflectionPrevious, positionsInPreviousTimeStep);
	GetDeflectionVectorFromPositions(deflectionNow, nodePositionsRelativeToRoot);
	NewmarkSolve(u2Out, newmarkMatrixR1CholeskyDecomposed, newmarkMatrixR2, newmarkMatrixR3, load, deflectionPrevious, deflectionNow);
}

void FemDeformation::CreateBeamSections()
{
#ifdef _DEBUG
	if (!beamSections.empty())
	{
		std::cout << "Valve already had BeamSections. Deleting the ones that are here now." <<std::endl; 
		beamSections.clear();
	}
#endif
	

	double currentNodePosX = 0;
	double length, leftWidth, rightWidth, leftThickness, rightThickness;
	bool bIsFixed;

	for (int nodeIndex = 0; nodeIndex < amountOfNodes - 1; nodeIndex++) //todo: check if this actually sets the fixed nodes properly.
	{
		// The 'fixed' sections are always considered at constant properties.
		if (nodeIndex < fixedNodes)
		{
			length = fixedLength / fixedNodes;
			leftWidth = rootWidth;
			rightWidth = rootWidth;
			leftThickness = rootThickness;
			rightThickness = rootThickness;
			bIsFixed = true;
		}
		else
		{
			bIsFixed = false;
			const double ratioCovered = currentNodePosX / freeLength;
			// It's continuous with the previous element, so re-use those values!
			leftWidth = beamSections.back().b[1];
			leftThickness = beamSections.back().h[1];
			// Handling how the free element is created based on the selected profile.
			switch (geometry.beamProfile)
			{
			case EBeamProfile::STRAIGHT_DOUBLE_TAPERED:
				rightWidth = rootWidth + (tipWidth - rootWidth)*ratioCovered;
				rightThickness = rootThickness + (tipThickness - rootThickness)*ratioCovered;
				break;
			default:
				throw std::logic_error("Determining beam properties is not implemented for this type of beam profile!");
			}
			// TODO: Right now, all nodes are exactly equidistant, and the same size. Make this variable with a spacing parameter.
			length = freeLength / freeNodes;
		}
		bool bCalculatePressureLoad = !bIsFixed; // if it is considered 'fixed', it means that one or two of the nodes it is connected to cannot move. Then don't calculate a pressure load, as this load will be distributed to the two nodes around it.
		BeamSection beamSection = BeamSection(length, {leftWidth, rightWidth}, {leftThickness, rightThickness}, material.density, material.youngsModulus, geometry.beamProfile, bIsFixed, nodeIndex ); //If you want material density and youngsModulus to be variable over the different sections, this is where you need to propogate the values to.
		beamSections.emplace_back(beamSection);

		// Advancing & resetting the values so it's easier to spot if there's a mistake.
		currentNodePosX += length;

		#ifdef _DEBUG
		leftWidth = 0;
		rightWidth = 0;
		leftThickness = 0;
		rightThickness = 0;
		#endif
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
			dampingMatrix(i,j) = empiricalParameters.rayleighDampingAlpha * globalMassMatrix(i,j) + empiricalParameters.rayleighDampingBeta * globalStiffnessMatrix(i,j);
		}
	}
	
	#ifdef _DEBUG
	assert(matrixOut.HasDiagonalGrainsOnly(4));
	#endif
}

void FemDeformation::AssembleNewmarkMatrix(TwoDimensionalArray& R1CholeskyOut, TwoDimensionalArray& R2Out, TwoDimensionalArray& R3Out, TwoDimensionalArray& KCholeskyOut, const double dt) const
{
	if (globalStiffnessMatrix.IsEmpty() || globalMassMatrix.IsEmpty() || dampingMatrix.IsEmpty())
		throw std::logic_error("Cannot assemble Newmark matrix, as the stiffness matrix or mass matrix has not yet been initialised.");

	R2Out.Resize(N_DOF, N_DOF);
	R3Out.Resize(N_DOF, N_DOF);

	// We will do a cholesky decomposition + transposition on R1, so separate it for now.
	TwoDimensionalArray R1PreProcessing(N_DOF,N_DOF, 0);
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

std::vector<int> FemDeformation::GetDOFVector() const
{
	int dofIndex = 0;
	std::vector<int> DOFVector(freeNodes, 0);
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
void FemDeformation::GetDeflectionVectorFromPositions(std::vector<double>& deflectionVectorOut, const std::vector<Position>& positions)
{
	deflectionVectorOut.resize(positions.size()*2);
	for (int i = 0; i < positions.size(); i++)
	{
		//TODO: Doublecheck if the deflection vectors are in fact x,y,x,y and not y,x,y,x
		deflectionVectorOut[2*i] = positions[i].x;
		deflectionVectorOut[2*i + 1] = positions[i].y;
	}
}

/* Given a symmetric positive definite matrix M (size NxN) (should be checked by user), this function computes the Cholesky decomposition of this matrix and fills a matrix L (that should be allocated and initialized by the user) so that A = L*LT (LT is the transpose matrix of L). The output L matrix is a superior triangular matrix. */
TwoDimensionalArray FemDeformation::CholeskyDecomposition(const TwoDimensionalArray& matrix, const std::vector<int>& DOFVector)
{
	#ifdef _DEBUG
	// Technically, the matrix m ust be positive-definite, which is a subset of symmetric matrices where all the local determinants are positive. However, testing that is relatively costly for a big matrix and its a big long algorithm, so it's left out- this is just to catch the obvious cases anyway!
	if (!matrix.IsDiagonallySymmetric())
		throw std::logic_error("An matrix that is not diagonally symmetric cannot be cholesky decomposed.");
	#endif
	
	// This is a variation on the Cholesky-Crout algorithm?
	//const int activeDegreesOfFreedom = freeNodes * N_DOF_PER_NODE; <- the original thing, but should work like this too.
	const int activeDegreesOfFreedom = matrix.nX;
	TwoDimensionalArray out(matrix.nX, matrix.nX, 0);
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
		out(i,i) = std::sqrt(out(i, i));

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
void FemDeformation::NewmarkSolve(std::vector<double>& u2Out, const TwoDimensionalArray& R1Cholesky, const TwoDimensionalArray& R2, const TwoDimensionalArray& R3, const std::vector<double>& load, const std::vector<double>& u0, const std::vector<double>& u1) const
{
	#ifdef _DEBUG
	assert(load.size() == u2Out.size());
	assert(u0.size() == u1.size());
	assert(load.size() == u1.size());
	assert(static_cast<int>(load.size()) == R1Cholesky.nX);
	assert(R1Cholesky.nX == R1Cholesky.nY);
	assert(R2.nX == R2.nY);
	assert(R3.nX == R3.nX);
	assert(R1Cholesky.nX == R2.nX);
	assert(R2.nX == R3.nX);
	#endif
	
	const std::vector<int> activeDof = GetDOFVector();
	double sum;
	std::vector<double> b = load;
	
	// Compute right-hand term of "A1*U(n+1) = F(n) + A2*U(n) + A3*U(n-1)"
	for (int i = 0; i < N_DOF; ++i)
	{
		for (int j = 0; j < N_DOF; ++j)
		{
			b[i] += R2.GetAt(i,j)*u1[j] + R3.GetAt(i,j)*u0[j];
		}
	}

	// First step : Forward substitution starting from first index
	for (int i = 0; i < n_active; ++i)
	{
		sum = 0.0;
		for (int j = 0; j < i; ++j)
		{
			sum += R1Cholesky.GetAt(i,j)*u2Out[activeDof[j]];
		}
		u2Out[activeDof[i]] = (b[activeDof[i]] - sum)/R1Cholesky.GetAt(i,i);
	}

	// Second step : Backward subsitution starting from last index
	for (int i = n_active - 1; i > -1; --i)
	{
		sum = 0.0;
		for (int j = i + 1; j < n_active; ++j)
		{
			sum += R1Cholesky.GetAt(j,i)*u2Out[activeDof[j]];
		}
		u2Out[activeDof[i]] = (u2Out[activeDof[i]] - sum)/R1Cholesky.GetAt(i,i);
	}

	#ifdef _DEBUG
	for (int i = 0; i < N_DOF; ++i)
	{
		if (std::isnan(u2Out[i]))
			throw std::logic_error("Encountered a NaN in newmark solving!");
	}
	#endif

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
	std::vector<int> DOFVector = GetDOFVector();


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
		deflectionVectorOut[DOFVector.at(i)] = (deflectionVectorOut[DOFVector.at(i)] - sum)/globalStiffnessMatrixCholeskyDecomposed.GetAt(i,i);
	}

	#ifdef _DEBUG
	for (const auto num : deflectionVectorOut)
		if(std::isnan(num))
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
Position FemDeformation::GetPositionOfBeamSection(const BeamSection &beamSection) const
{
	return (nodePositionsRelativeToRoot.at(beamSection.leftNodeIndex) + nodePositionsRelativeToRoot.at(beamSection.rightNodeIndex)) * 0.5;
}