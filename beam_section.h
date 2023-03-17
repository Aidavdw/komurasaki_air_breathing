#include <array>

struct BeamSection
{
	struct BeamSection(const int id, const double posX, const double length, const double width[2], const double thickness[2], const double density);

	int id;

	double posX;					// Position of the bottom left corner of the beam section relative to the root.
	double b[2];					// The width of the beam section's root part. Since this is a 2d case, this is in the circumferential direction (into the paper, so-to-speak).
	double h[2];					// Thickness of the beam section. technically the 'height', therefore the odd name choice. Might change later.
	double length;

	double density;


	double crossSectionalArea[2];	// cross-sectional area of the beam section, orthogonal to the bending plane.
	double areaMomentOfInertia[2];	// the area moment of inertia I, in m^4.
	double massMatrix[4][4];

};