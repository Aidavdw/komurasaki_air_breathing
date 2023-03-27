#pragma once

// Represents the properties of a beam section in a fem simulation. Note that it does not say anything about its deformation state- these are all just constant properties.
struct BeamSection
{
	struct BeamSection(const double length, const double width[2], const double thickness[2], const double density, const double youngsModulus);

	double b[2];					// The width of the beam section's root part. Since this is a 2d case, this is in the circumferential direction (into the paper, so-to-speak).
	double h[2];					// Thickness of the beam section. technically the 'height', therefore the odd name choice. Might change later.
	double length;

	double density;
	double youngsModulus;			// The Young's modulus E expressed in Pa(scals)


	double crossSectionalArea[2];	// cross-sectional area of the beam section, orthogonal to the bending plane.
	double areaMomentOfInertia[2];	// the area moment of inertia I, in m^4.
	double massMatrix[4][4];
	double stiffnessMatrix[4][4];

private:
	void PopulateMassMatrix();
	void PopulateStiffnessMatrix();

};