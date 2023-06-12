﻿#include "main_library.h"
#include "sim_case.h"
#include "extern/pybind11/include/pybind11/pybind11.h"

namespace py = pybind11;

PYBIND11_MODULE(komurasakiairbreathing, m)
{
    // Simcase
    py::class_<SimCase>(m, "SimCase")
        .def(py::init<int, int>()) // constructor
        .def("AddDomain", &SimCase::AddDomain)
        .def("ConnectBoundariesByName", &SimCase::ConnectBoundariesByName)
        .def("ConnectBoundariesById", &SimCase::ConnectBoundariesById);

    // Todo: add these jetsers
    // SolverSettings
    // AmbientConditions
    // ChapmanJougetInitialConditionParameters
    // RuntimeParameters

    // Position
    py::class_<Position>(m,"Position")
        .def(py::init<int, int>()) // constructor
        .def_readonly("x", &Position::x)
        .def_readonly("y", &Position::y)
        .def_readonly("up_direction", &Position::upDirection);
    

    // MeshSpacing
    py::class_<MeshSpacing>(m,"MeshSpacing")
        .def(py::init<EMeshSpacingType, double, int, double, double>()) // constructor
        .def_readonly("left", &MeshSpacing::left)
        .def_readonly("right", &MeshSpacing::right)
        .def_readonly("length", &MeshSpacing::length)
        .def_readonly("amount_of_elements", &MeshSpacing::amountOfElements);

    // EInitialisationMethod
    py::enum_<EInitialisationMethod>(m, "InitialisationMethod")
        .value("zero",EInitialisationMethod::ZERO)
        .value("ambient_conditions",EInitialisationMethod::AMBIENT_CONDITIONS)
        .value("from_input_data",EInitialisationMethod::FROM_INPUT_DATA)
        .value("from_chapman_jouget_solution",EInitialisationMethod::FROM_CHAPMAN_JOUGET_SOLUTION)
        .export_values();

    // EFace
    py::enum_<EFace>(m, "Face")
        .value("left",EFace::LEFT)
        .value("right",EFace::RIGHT)
        .value("top",EFace::TOP)
        .value("bottom",EFace::BOTTOM)
        .export_values();

    // EBoundaryCondition
    py::enum_<EBoundaryCondition>(m, "BoundaryCondition")
        .value("not_set",EBoundaryCondition::NOT_SET)
        .value("slip",EBoundaryCondition::SLIP)
        .value("no_slip",EBoundaryCondition::NO_SLIP)
        .value("connected",EBoundaryCondition::CONNECTED)
        .value("supersonic_inlet",EBoundaryCondition::SUPERSONIC_INLET)
        .value("supersonic_outlet",EBoundaryCondition::SUPERSONIC_OUTLET)
        .export_values();

    // Mesh Spacing Type
    py::enum_<EMeshSpacingType>(m, "MeshSpacingType")
        .value("constant",EMeshSpacingType::CONSTANT)
        .value("linear",EMeshSpacingType::LINEAR)
        .value("parabolic",EMeshSpacingType::PARABOLIC)
        .value("exponential",EMeshSpacingType::EXPONENTIAL)
        .export_values();
    

    // define all standalone functions
    m.def("DoSimulation", &DoSimulation);
}
