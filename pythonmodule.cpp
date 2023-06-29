#include "main_library.h"
#include "reed_valve.h"
#include "sim_case.h"
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "pythoninterface/conversion_functions.h"
#include "extern/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(komurasakiairbreathing, m)
{
    // Simcase
    py::class_<SimCase>(m, "SimCase")
        .def(py::init<>()) // default constructor
        .def("AddDomain", &SimCase::AddDomain)
        .def("ConnectBoundariesByName", &SimCase::ConnectBoundariesByName)
        .def("ConnectBoundariesById", &SimCase::ConnectBoundariesById)
        .def("GetDomainByID", &SimCase::GetDomainByID)
        .def("GetDomainByName", &SimCase::GetDomainByName)
        .def("AddRecord", &SimCase::AddRecord)
        .def_readonly("domains", &SimCase::domains)
        .def_readonly("two_dimensional_array_records", &SimCase::twoDimensionalArrayRecords)
        .def_readwrite("dt", &SimCase::dt)
        .def_readwrite("simulation_duration", &SimCase::simulationDuration);
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

    py::class_<ReedValve>(m, "ReedValve")
        .def_readonly("pos", &ReedValve::pos)
        .def_readonly("hole_end_position_in_domain", &ReedValve::holeEndPositionInDomain);

    py::class_<TwoDimensionalArray>(m, "TwoDimensionalArray")
        .def("GetAt", &TwoDimensionalArray::GetAtPythonProxy)
        .def_readonly("size_x", &TwoDimensionalArray::nX)
        .def_readonly("size_y", &TwoDimensionalArray::nY);

    py::class_<MaterialProperties>(m, "MaterialProperties")
        .def(py::init<>()) // default constructor
        .def_readwrite("youngs_modulus", &MaterialProperties::youngsModulus)
        .def_readwrite("density", &MaterialProperties::density);

    py::class_<ReedValveEmpiricalParameters>(m, "ReedValveEmpiricalParameters")
        .def(py::init<>()) // default constructor
        .def_readwrite("natural_frequency", &ReedValveEmpiricalParameters::naturalFrequency)
        .def_readwrite("rayleigh_damping_alpha", &ReedValveEmpiricalParameters::rayleighDampingAlpha)
        .def_readwrite("rayleigh_damping_beta", &ReedValveEmpiricalParameters::rayleighDampingBeta)
        .def_readwrite("aerodynamic_damping_c1", &ReedValveEmpiricalParameters::dampingC1)
        .def_readwrite("aerodynamic_damping_c2", &ReedValveEmpiricalParameters::dampingC2)
        .def_readwrite("aerodynamic_damping_c3", &ReedValveEmpiricalParameters::dampingC3)
        .def_readwrite("hole_factor", &ReedValveEmpiricalParameters::holeFactor);

    py::class_<ReedValveGeometry>(m, "ReedValveGeometry")
        .def(py::init<>()) // default constructor
        .def_readwrite("free_length", &ReedValveGeometry::freeLength)
        .def_readwrite("root_thickness", &ReedValveGeometry::rootThickness)
        .def_readwrite("tip_thickness", &ReedValveGeometry::tipThickness)
        .def_readwrite("root_width", &ReedValveGeometry::rootWidth)
        .def_readwrite("tip_width", &ReedValveGeometry::tipWidth);

    py::class_<SolverSettings>(m, "SolverSettings")
        .def(py::init<>()) // default constructor
        .def_readwrite("ausm_switch_bias", &SolverSettings::AUSMSwitchBias)
        .def_readwrite("muscl_bias", &SolverSettings::MUSCLBias)
        .def_readwrite("entropy_fix", &SolverSettings::entropyFix)
        .def_readwrite("runge_kutta_order", &SolverSettings::rungeKuttaOrder)
        .def_readwrite("amount_of_ghost_cells", &SolverSettings::nGhost)
        .def_readwrite("flux_limiter_type", &SolverSettings::nGhost);

    py::enum_<EFluxLimiterType>(m, "FluxLimiterType")
        .value("none",EFluxLimiterType::NONE)
        .value("min_mod",EFluxLimiterType::MIN_MOD)
        .value("super_bee",EFluxLimiterType::SUPER_BEE)
        .value("van_albada_one",EFluxLimiterType::VAN_ALBADA_ONE)
        .value("van_albada_two",EFluxLimiterType::VAN_ALBADA_TWO)
        .export_values();
    
    py::class_<AmbientConditions>(m, "AmbientConditions")
        .def(py::init<>()) // default constructor
        .def_readwrite("mach",&AmbientConditions::mach)
        .def_readwrite("temperature",&AmbientConditions::temperature)
        .def_readwrite("static_pressure",&AmbientConditions::staticPressure)
        .def_readwrite("free_flow_x_velocity",&AmbientConditions::u)
        .def_readwrite("free_flow_y_velocity",&AmbientConditions::v);

    py::class_<ChapmanJougetInitialConditionParameters>(m, "ChapmanJougetInitialConditionParameters")
        .def(py::init<>()) // default constructor
        .def_readwrite("beam_power",&ChapmanJougetInitialConditionParameters::beamPower)
        .def_readwrite("energy_absorption_coefficient",&ChapmanJougetInitialConditionParameters::energyAbsorptionCoefficient)
        .def_readwrite("specific_heat_ratio",&ChapmanJougetInitialConditionParameters::gamma)
        .def_readwrite("ideal_gas_constant",&ChapmanJougetInitialConditionParameters::idealGasConstant);

    py::class_<RuntimeParameters>(m, "RuntimeParameters")
        .def(py::init<>()) // default constructor
        .def_readwrite("time_steps_between_data_export", &RuntimeParameters::numberOfTimeStepsBetweenDataExport);

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
        .def(py::init<>()) // constructor
        .def_readwrite("type", &MeshSpacing::spacingType)
        .def_readwrite("left", &MeshSpacing::left)
        .def_readwrite("right", &MeshSpacing::right)
        .def_readwrite("amount_of_elements", &MeshSpacing::amountOfElements);

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

    py::enum_<EDirection>(m, "Direction")
        .value("right", EDirection::RIGHT)
        .value("up", EDirection::UP)
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
        .value("straight_double_tapered",EBeamProfile::STRAIGHT_DOUBLE_TAPERED)
        .export_values();
    

    // define all standalone functions
    m.def("DoSimulation", &DoSimulation);
    m.def("TwoDimensionalArrayToNumpyArray", &TwoDimensionalArrayToNumpyArray);
    m.def("FillTwoDimensionalArrayFromNumpy", &FillTwoDimensionalArrayFromNumpy);
}
