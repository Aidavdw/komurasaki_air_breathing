#include "valve.h"
#include "fem_deformation.h"
#include <vector>


// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public Valve, public FemDeformation
{
public:
	ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections);

	std::vector<int[2]> sourceCellIndices;

	void OnRegister() override;

private:
	// 
	void SetSourceCellIndices(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections);
	void SetPressureReadingCellIndices(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections);
};