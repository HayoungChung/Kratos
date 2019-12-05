// System includes

#if defined(KRATOS_PYTHON)
// External includes


// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "analytic_sensitivity_application.h"
#include "./custom_elements/small_displacement_ersatz_element.hpp"
#include "./custom_strategies/structure_adjoint_sensitivity_strategy.h"
#include "./custom_strategies/custom_schemes/residualbased_incrementalupdate_static_ersatz_scheme.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosAnalyticSensitivityApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosAnalyticSensitivityApplication,
            KratosAnalyticSensitivityApplication::Pointer,
            KratosApplication >(m,"KratosAnalyticSensitivityApplication")
            .def(py::init<>())
            ;

    // AddCustomStrategiesToPython(m);
    // AddCustomProcessesToPython(m);
    // AddCustomUtilitiesToPython(m);
    // AddCustomConstitutiveLawsToPython(m);
    // AddCustomResponseFunctionUtilitiesToPython(m);

    // py::class_<Variable<ShellCrossSection::Pointer>,VariableData >(m,"ShellCrossSectionVariable");

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, RHO_MIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, X_PHYS ) // material density
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, X_PHYS_OLD ) 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DFDX ) // sensitivity in an element
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DFDX_intg ) // sensitivity in an element
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, SOLID_VOID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LOCAL_STRAIN_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, COMPLIANCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, COMPLIANCE_intg )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, VOLUME )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, P_NORM_STRESS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, UF )

    // py::class_< SmallDisplacementErsatzElement, SmallDisplacementErsatzElement::Pointer, Element >  // without , Element>, the methods exposed to Element.h (mother class) is not usable in this element
    // (m, "SmallDisplacementErsatzElement")
    //     .def(py::init<>())
    //     .def("GetValueOnIntegrationPoints", &SmallDisplacementErsatzElement::GetValueOnIntegrationPoints)
    //     .def("Calculate", &SmallDisplacementErsatzElement::Calculate)
    //     .def("CalculateOnIntegrationPoints", &SmallDisplacementErsatzElement::CalculateOnIntegrationPoints)
    //     ;

	// py::class_< ResidualBasedIncrementalUpdateStaticErsatzchemeType, ResidualBasedIncrementalUpdateStaticErsatzchemeType::Pointer, ResidualBasedIncrementalUpdateStaticScheme >
	// (m, "ResidualBasedIncrementalUpdateStaticErsatzcheme")
    // .def(py::init<>())
	// .def("Initialize", &ResidualBasedIncrementalUpdateStaticErsatzcheme<SparseSpaceType, LocalSpaceType>::Initialize)
	// ;}

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef StructureAdjointSensitivityStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > StructureAdjointSensitivityStrategyType;
	typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;

    py::class_< StructureAdjointSensitivityStrategyType,StructureAdjointSensitivityStrategyType::Pointer,  BaseSolvingStrategyType >
	(m, "StructureAdjointSensitivityStrategy")
    .def(py::init<ModelPart&, LinearSolverType::Pointer, int>())
	.def("ComputeComplianceSensitivities_elem",&StructureAdjointSensitivityStrategyType::ComputeComplianceSensitivities_elem)
	.def("ComputeVolumeFractionSensitivities_elem",&StructureAdjointSensitivityStrategyType::ComputeVolumeFractionSensitivities_elem)
    .def("ComputeComplianceSensitivities_gpt",&StructureAdjointSensitivityStrategyType::ComputeComplianceSensitivities_gpt)
	.def("ComputeVolumeFractionSensitivities_gpt",&StructureAdjointSensitivityStrategyType::ComputeVolumeFractionSensitivities_gpt)
    .def("ComputeGaussPointCoordinates",&StructureAdjointSensitivityStrategyType::ComputeGaussPointCoordinates)
	;

	typedef ResidualBasedIncrementalUpdateStaticErsatzScheme< SparseSpaceType, LocalSpaceType > ResidualBasedIncrementalUpdateStaticErsatzSchemeType;
	typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    py::class_< ResidualBasedIncrementalUpdateStaticErsatzSchemeType, ResidualBasedIncrementalUpdateStaticErsatzSchemeType::Pointer, BaseSchemeType >
	(m, "ResidualBasedIncrementalUpdateStaticErsatzScheme")
    .def(py::init<>())
	.def("Initialize", &ResidualBasedIncrementalUpdateStaticErsatzScheme<SparseSpaceType, LocalSpaceType>::Initialize)
	;
}
}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
