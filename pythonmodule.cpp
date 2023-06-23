#include "main_library.h"
#include "sim_case.h"
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "pythoninterface/conversion_functions.h"
#include "extern/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(komurasakiairbreathing, m)
{
    // Simcase
    py::class_<SimCase>(m, "SimCase")
        .def(py::init<double, double>()) // constructor
        .def("AddDomain", &SimCase::AddDomain)
        .def("ConnectBoundariesByName", &SimCase::ConnectBoundariesByName)
        .def("ConnectBoundariesById", &SimCase::ConnectBoundariesById)
        .def("GetDomainByID", &SimCase::GetDomainByID)
        .def("GetDomainByName", &SimCase::GetDomainByName)
        .def("AddRecord", &SimCase::AddRecord)
        .def_readonly("domains", &SimCase::domains)
        .def_readonly("two_dimensional_array_records", &SimCase::twoDimensionalArrayRecords);
        //.def_readonly("valve_interfaces", &SimCase::valves);
    

    py::class_<Domain>(m, "Domain")
        .def("SetBoundaryType", &Domain::SetBoundaryType)
        .def_readonly("density", &Domain::rho)
        .def_readonly("velocity_x", &Domain::u)
        .def_readonly("velocity_r", &Domain::v)
        .def_readonly("pressure", &Domain::p)
        .def_readonly("internal_energy", &Domain::E)
        .def_readonly("temperature", &Domain::T)
        .def_readonly("enthalpy", &Domain::H);

    py::class_<TwoDimensionalArrayRecord>(m,"TwoDimensionalArrayRecord")
        .def_readonly("records", &TwoDimensionalArrayRecord::records)
        .def("AsNumpyArray", &TwoDimensionalArrayRecord::AsNumpyArray);
        

    py::class_<FieldQuantity>(m, "FieldQuantity")
        .def_readonly("current_time_step", &FieldQuantity::currentTimeStep)
        .def_readonly("rungeKuttaBuffer", &FieldQuantity::currentTimeStep);

    py::class_<TwoDimensionalArray>(m, "TwoDimensionalArray")
        .def("GetAt", &TwoDimensionalArray::GetAtPythonProxy)
        .def_readonly("size_x", &TwoDimensionalArray::nX)
        .def_readonly("size_y", &TwoDimensionalArray::nY);

    py::class_<MaterialProperties>(m, "MaterialProperties")
        .def(py::init<MaterialProperties>()) // constructor
        .def_readwrite("youngs_modulus", &MaterialProperties::youngsModulus)
        .def_readwrite("density", &MaterialProperties::density);

    py::class_<ReedValveEmpiricalParameters>(m, "ReedValveEmpiricalParameters")
        .def(py::init<ReedValveEmpiricalParameters>()) // constructor
        .def_readwrite("natural_frequency", &ReedValveEmpiricalParameters::naturalFrequency)
        .def_readwrite("rayleigh_damping_alpha", &ReedValveEmpiricalParameters::rayleighDampingAlpha)
        .def_readwrite("rayleigh_damping_beta", &ReedValveEmpiricalParameters::rayleighDampingBeta)
        .def_readwrite("aerodynamic_damping_c1", &ReedValveEmpiricalParameters::dampingC1)
        .def_readwrite("aerodynamic_damping_c2", &ReedValveEmpiricalParameters::dampingC2)
        .def_readwrite("aerodynamic_damping_c3", &ReedValveEmpiricalParameters::dampingC3)
        .def_readwrite("hole_factor", &ReedValveEmpiricalParameters::holeFactor);

    py::class_<ReedValveGeometry>(m, "ReedValveGeometry")
        .def(py::init<ReedValveGeometry>()) // constructor
        .def_readwrite("free_length", &ReedValveGeometry::freeLength)
        .def_readwrite("root_thickness", &ReedValveGeometry::rootThickness)
        .def_readwrite("tip_thickness", &ReedValveGeometry::tipThickness)
        .def_readwrite("root_width", &ReedValveGeometry::rootWidth)
        .def_readwrite("tip_width", &ReedValveGeometry::tipWidth);
    

    // Todo: add these jetsers
    // SolverSettings
    // AmbientConditions
    // ChapmanJougetInitialConditionParameters
    // RuntimeParameters

    py::class_<Boundary>(m,"Boundary")
        .def_readonly("boundary_type", &Boundary::boundaryType);

    // Position
    py::class_<Position>(m,"Position")
        .def(py::init<double, double>()) // constructor
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
        .value("from_input_data",EInitialisationMethod::FROM_INPUT_RHO_P_U_V)
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

    py::enum_<EBeamProfile>(m, "BeamProfile")
    .value("straigh_double_tapered",EBeamProfile::STRAIGHT_DOUBLE_TAPERED)
    .export_values();
    

    // define all standalone functions
    m.def("DoSimulation", &DoSimulation);
    m.def("TwoDimensionalArrayToNumpyArray", &TwoDimensionalArrayToNumpyArray);
    m.def("FillTwoDimensionalArrayFromNumpy", &FillTwoDimensionalArrayFromNumpy);
}
