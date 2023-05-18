#pragma once
#include <utility>

// Describes what the beam looks like; is it straight? is it double tapered?
enum class EBeamProfile
{
	STRAIGHT_DOUBLE_TAPERED // Has variable width (b) and variable thickness/height (h). Is exactly straight.
};

// Represents the properties of a beam section in a fem simulation. Note that it does not say anything about its deformation state- these are all just constant properties.
class BeamSection
{
public:
	BeamSection(const double length, const std::pair<double,double> width, const std::pair<double,double> thickness, const double density, const double youngsModulus, const EBeamProfile beamProfile, const bool bIsFixed, const int leftNodeIndex);

	int leftNodeIndex;
	int rightNodeIndex;

	EBeamProfile beamProfile;
	bool bHasPressureLoad;
	

	double b[2] = {0,0};					// The width of the beam section's root part. Since this is a 2d case, this is in the circumferential direction (into the paper, so-to-speak).
	double h[2] = {0,0};					// Thickness of the beam section. technically the 'height', therefore the odd name choice. Might change later.
	double length;
	double topOrBottomSurfaceArea; // the area on top of (or on the bottom of) the beam section. Is used to determine the loads in combination with the given pressure difference.

	double density;
	double youngsModulus;			// The Young's modulus E expressed in Pa(scals)

	double crossSectionalArea[2]{};	// cross-sectional area of the beam section, orthogonal to the bending plane.
	double areaMomentOfInertia[2]{};	// the area moment of inertia I, in m^4.
	double massMatrix[4][4]{};
	double stiffnessMatrix[4][4]{};

private:
	void PopulateMassMatrix();
	void PopulateStiffnessMatrix();

};