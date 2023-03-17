#include <vector>
#include "2dArray.h"
#include "beam_section.h"

// Describes what the beam looks like; is it straight? is it double tapered?
enum EBeamProfile
{
	STRAIGHTDOUBLETAPERED // Has variable width (b) and variable thickness/height (h). Is exactly straight.
};

// Abstract class that implements data and methods for simple 1d fem deformation calculations
class FemDeformation
{

	class FemDeformation(const int amountOfFreeSections, const int amountOfFixedNodes, const EBeamProfile beamProfile);

	EBeamProfile beamProfile;
	std::vector<BeamSection> beamSections;		// The individual FEM segments in this valve. Note that this is one less than there are nodes!
	int fixedNodes;								// The amount of sections in the beam that are considered 'fixed'
	int amountOfNodes;							// The total amount of nodes that this beam is modeled with. This means fixed, and free nodes.
	int N_DOF;									// The amoutn of degrees of freedom for the FEM system.

	double freeLength;							// Length of the part that can move freely
	double fixedLength;							// Length of the part that is fixed in place.

	double rootWidth;
	double tipWidth;
	double rootThickness;
	double tipThickness;


private:
	void CreateBeamSections();
	void PopulateGlobalMassMatrix(TwoDimensionalArray& matrixOut);
	void PopulateGlobalStiffnessMatrix(TwoDimensionalArray& matrixOut);
};