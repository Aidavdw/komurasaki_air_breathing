
struct BeamSection
{
	struct BeamSection(const double posX, const double width, const double thickness);

	double posX;				// Position of the bottom left corner of the beam section relative to the root.
	double b;					// The width of the beam section. Since this is a 2d case, this is in the circumferential direction (into the paper, so-to-speak).
	double h;					// Thickness of the beam section. technically the 'height', therefore the odd name choice. Might change later.

	double crossSectionalArea;	// cross-sectional area of the beam section, orthogonal to the bending plane.
	double areaMomentOfInertia;	// the area moment of inertia I, in m^4.

};