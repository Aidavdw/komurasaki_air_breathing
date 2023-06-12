#pragma once

#include <../extern/pybind11/include/pybind11/pybind11.h>
#include <../extern/pybind11/include/pybind11/numpy.h>

#include "2dArray.h"
namespace py = pybind11;

py::array_t<double> TwoDimensionalArrayToNumpyArray(const TwoDimensionalArray& src);