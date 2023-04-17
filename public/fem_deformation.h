#pragma once
#include <vector>
#include "2dArray.h"
#include "beam_section.h"
#include "pos2d.h"



// A simple beam in bending, modeled as an FEM system.
class FemDeformation
{
public:
	FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const EBeamProfile beamProfile, const double freeLength, const double fixedLength, const double dt);

	FemDeformation();
	
	std::vector<BeamSection> beamSections;		// Data on the individual FEM segments in this valve. Note that this is one less than there are nodes!
	int fixedNodes;								// The amount of nodes (=sections+1) in the beam that are considered 'fixed', unable to deform.
	int freeNodes;								//  The amount of nodes (=sections+1) in the beam that are considered able to deform
	int amountOfNodes;							// The total amount of nodes that this beam is modeled with. This means fixed, and free nodes.
	int N_DOF;									// The amoutn of degrees of freedom for the FEM system.

	double dt;									// Time step used to assemble the newmark matrices.
	
	std::vector<Position> nodePositionsRelativeToRoot;	// The current positions of the nodes, in a local coordinate space.
	std::vector<Position> positionsInPreviousTimeStep;
	
protected:
	EBeamProfile beamProfile_;
	double freeLength;							// Length of the part that can move freely
	double fixedLength;							// Length of the part that is fixed in place.

	double rayleighDampingAlpha;				 // Alpha coef. for Rayleigh damping (alpha*M + beta*K)
	double rayleighDampingBeta;				 // Beta coef. for Rayleigh damping (alpha*M + beta*K)
	
	double rootWidth;
	double tipWidth;
	double rootThickness;
	double tipThickness;

	TwoDimensionalArray globalMassMatrix;		// The global mass matrix for this FEM beam, with the element matrices combined.
	TwoDimensionalArray globalStiffnessMatrix;	// The global mass stiffness for this FEM beam, with the element matrices combined.
	TwoDimensionalArray globalStiffnessMatrixCholeskyDecomposed;
	TwoDimensionalArray dampingMatrix;			// The global mass stiffness for this FEM beam, with the element matrices combined.
	
	// Newmark matrices
	TwoDimensionalArray newmarkMatrixR1CholeskyDecomposed;	// The first newmark matrix, but cholesky decomposed. 
	TwoDimensionalArray newmarkMatrixR2;					// The second newmark matrix.
	TwoDimensionalArray newmarkMatrixR3;					// The third newmark matrix.

public:

	// Updates the position of the fem elements based on the current pressure field.
	// todo: figure out what u1 and u2 actually are.
	void UpdatePositions(const std::vector<double>& u1, const std::vector<double>& u2);
	
	// Solves the system of equations for the stiffness of the cholesky-decomposed stiffness matrix. Is only used once to set up the problem.
	void SolveCholeskySystem(std::vector<double>& deflectionVectorOut, const std::vector<double>& load) const;

	// Every node (bar the first one and the last one) has two beam sectioned connected to it. This gets a reference to those two beam sections. The optional argument says whether to return the left (false) value or the right (true) value
	BeamSection* BeamSectionsConnectedToNode(const int nodeIndex, const bool bRight);

	



protected:
	// Creates geometry for the nodes and the edges between them. Populates the beamSections and nodePositionsRelativeToRoot.
	void CreateBeamSections();
	void AssembleGlobalMassMatrix(TwoDimensionalArray& matrixOut) const;
	void AssembleGlobalStiffnessMatrix(TwoDimensionalArray& matrixOut) const;

	std::vector<double> GetDOFVector() const;

	// Build Damping matrix based on Rayleigh's damping model. ALPHA and BETA factors (respectively of M and K) should be specified by the user. These coefficients can be known from experiment.
	void AssembleDampingMatrix(TwoDimensionalArray& matrixOut);

	// Computes matrices needed for resolution according to Newmark's scheme
	void AssembleNewmarkMatrix(TwoDimensionalArray& R1CholeskyOut, TwoDimensionalArray& R2Out, TwoDimensionalArray& R3Out, TwoDimensionalArray& KCholeskyOut, const double dt) const;

	static TwoDimensionalArray CholeskyDecomposition(const TwoDimensionalArray& matrix, const std::vector<double>& DOFVector);


};