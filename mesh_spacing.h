
struct Domain;


enum EMeshSpacingType
{
	CONSTANT, // Ignores everything, just assumes that all the cells are equally sized and distributed over the domain.
	LINEAR, // Scales the cells linearly from the left to the right. Ignores the centre parameter.
	PARABOLIC, // Sets cell sizes at the edges and at the center, and fits a parabolic function between them.
	EXPONENTIAL // Exponentially scales between the sides. Ignores the centre parameter.

};

// Contains information and calculations for the internal distribution of cells in a domain (mesh). Be sure to supply an adequate amount of spacing parameters for the desired mesh spacing type, as to not overconstrain the problem.
struct MeshSpacing
{
	// Could be implemented with polymorphism as well, but that would mean more lookups in the vtable based on virtual functions. Therefore, just state based on an enum.

	// Needs to have an empty constructor because of passing by val
	MeshSpacing()
	{};

	
	MeshSpacing(const EMeshSpacingType meshSpacingType, const double length, const int amountofElements, const double resolution_left, const double resolution_right);

	EMeshSpacingType spacingType = EMeshSpacingType::CONSTANT;
	
	double left = 0;
	double right = 0;

	double length = 0;
	int amountOfElements = 1;

	double GetCellWidth(const int i);

private:
	void FitSpacingToParameters();
	// Function describing how close to fitting the current configuration is
	double SpacingObjectiveFunction(std::vector<double>& funcLoc);
};




