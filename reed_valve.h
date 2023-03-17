#include "valve.h"
#include "fem_deformation.h"
#include <vector>


// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public Valve, public FemDeformation
{
	ReedValve(const int amountOfFreeSections, const int amountOfFixedNodes);

	void OnRegister() override;

};