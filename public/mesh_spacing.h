#pragma once
#include <vector>
#include "AuxFunctions.h"

// forward declarations
class Domain;

// For a domain, describes how the distribution of cell sizes changes. A domain has two of these, 1 in each axis.
enum class EMeshSpacingType
{
	CONSTANT, // Ignores everything, just assumes that all the cells are equally sized and distributed over the domain.
	LINEAR, // Scales the cells linearly from the left to the right. Ignores the centre parameter.
	PARABOLIC, // CURRENTLY UNIMPLEMENTED. Sets cell sizes at the edges and at the center, and fits a parabolic function between them.
	EXPONENTIAL // CURRENTLY UNIMPLEMENTED. Exponentially scales between the sides. Ignores the centre parameter.

};

// todo: remove double dependecy on 'length'. Make a proxy MeshSpacingWithoutSize, which gets constructed using the dimensions of the domain itself.

// Contains information and calculations for the internal distribution of cells in a domain (mesh). Be sure to supply an adequate amount of spacing parameters for the desired mesh spacing type, as to not overconstrain the problem.
struct MeshSpacing
{
	// Could be implemented with polymorphism as well, but that would mean more lookups in the vtable based on virtual functions. Therefore, just state based on an enum.

	// Needs to have an empty constructor because of passing by val
	MeshSpacing() = default;

	
	MeshSpacing(const EMeshSpacingType meshSpacingType, const double length, const int amountOfElements, const double resolutionLeft, const double resolutionRight);

	EMeshSpacingType spacingType = EMeshSpacingType::CONSTANT;
	
	double left = 0;
	double right = 0;

	double length = 0;
	int amountOfElements = 0;

	double GetCellWidth(const int i) const;

	bool operator== (const MeshSpacing& other) const
	{
		return (IsCloseToZero(left - other.left) && IsCloseToZero(right - other.right) && IsCloseToZero(length - other.length) && amountOfElements == other.amountOfElements);
	}

private:
	void FitSpacingToParameters();
	// Function describing how close to fitting the current configuration is
	double SpacingObjectiveFunction(std::vector<double>& funcLoc) const;
	
};




