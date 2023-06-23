#pragma once


// A small struct containing all the quantities in one cell. Differently organised from that in the domain, but this is a good way to transfer some data more easily.
struct CellValues
{
	CellValues() :
		density(0),
		u(0),
		v(0),
		p(0),
		e(0),
		t(0),
		h(0)
	{}

	double density; // Density
	double u; // velocity x-component
	double v; // velocity y-component
	double p; // pressure
	double e; // Internal energy?
	double t; // Temperature
	double h; // enthalpy ?
};