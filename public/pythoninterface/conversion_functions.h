#pragma once
/* This file contains functionality to convert from internal data representations to python representations for types that cannot be handled by PyBind11 directly. To keep code interdependency to a minimum, they are gathered here, and are not placed in their class definitions.
 */

#include <../extern/pybind11/include/pybind11/pybind11.h>
#include <../extern/pybind11/include/pybind11/numpy.h>

#include "2dArray.h"
#include "domain.h"
namespace py = pybind11;


py::array_t<double> TwoDimensionalArrayToNumpyArray(const TwoDimensionalArray& src); // Returns a copy of a two dimensional array as a python numpy array.

void SetDomainInitialValueFromNumpyArrays(Domain& domain, const py::array_t<double>& rho, const py::array_t<double>& p, const py::array_t<double>& u, const py::array_t<double>& v);

void FillTwoDimensionalArrayFromNumpy(TwoDimensionalArray& out, const py::array_t<double>& pyArray); // Sets a TwoDimensionalArray by reference based on a python array input