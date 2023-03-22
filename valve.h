#include "domain_enums.h"
#include <vector>

struct Domain;


//Abstract class that represents a valve that lets air from one domain into another.
class Valve
{
public:
	Valve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, double size);


	double size;
	EBoundaryLocation boundary;			// which side of the domain the valve is attached to
	double positionAlongBoundary;
	std::vector<double[2]> sourceIndices;
	std::vector<double[2]> pressureMeasuringIndices;

private:

	Domain* intoDomain;					// The domain that the valve creates a source into.

public:
	// Called to all valves upon being registered by the SimCase
	virtual void OnRegister();
	int GetSourceTermCount() const;

private:



};