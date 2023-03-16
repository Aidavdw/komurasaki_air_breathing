
// Abstract class that implements data and methods for simple 1d fem deformation calculations
class FemDeformation
{

	class FemDeformation(const int amountOfElementsToSplitBeamInto);

	int N_FEM;		// The amount of pieces that the beam is split up in
	int N_NODE;		// The total amount of connected nodes (points that can move, connected by elements) this beam has.
	int N_DOF;		// The total amount of degrees of freedom of this beam.
};