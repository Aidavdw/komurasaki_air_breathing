#include "finite_element.h"
#include <cmath>

BeamSection::BeamSection(const double posX, const double width, const double thickness) :
	posX(posX),
	b(width),
	h(thickness)
{
	crossSectionalArea = h * b;
	areaMomentOfInertia = b * pow(h, 3) / 12.0;
}