// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Hayoung Chung (Drived from topology optimization application)
//
// ==============================================================================

#if !defined(KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED)
#define  KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "../custom_elements/small_displacement_ersatz_element.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Application includes
#include "./analytic_sensitivity_application.h"


namespace Kratos {

///@addtogroup AnalyticSensitivityApplication
///@{

///@name Kratos Classes
///@{

/// Solution strategy to calculate the sensitivities.
/// Derives from the previously defined Solving Strategy

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class StructureAdjointSensitivityStrategy: public SolvingStrategy<TSparseSpace, TDenseSpace,TLinearSolver>
{
public:

	///@name Type Definitions
	///@{

	KRATOS_CLASS_POINTER_DEFINITION(StructureAdjointSensitivityStrategy);

	typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
	typedef typename Scheme<TSparseSpace,TDenseSpace>::Pointer SchemePointerType;
	typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderAndSolverPointerType;

	///@}
	///@name Life Cycle
	///@{

	StructureAdjointSensitivityStrategy( ModelPart& rStructureModelPart,
			typename TLinearSolver::Pointer pNewLinearSolver,
			const int dimension = 3)
	: BaseType(rStructureModelPart),
	  mr_structure_model_part(rStructureModelPart),
	  m_dimension(dimension)
	{}


	virtual ~StructureAdjointSensitivityStrategy()
	{}

	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- COMPUTE SENSITIVITIES  ------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// Computes compliance sensitivities from the adjoint solution (the problem is self-adjoint)
	void ComputeComplianceSensitivities_elem()
	{
		KRATOS_TRY;

		double Out = 0.0; // dummy variables

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			// Checking Calculate() - ok
			element_i->Calculate(COMPLIANCE, Out, ConstProcessInfo); 
			// std::cout << "Out = " << Out << d::endl;
			// std::cout << "GetValue = " << element_i->GetValue(COMPLIANCE) << std::endl;
			// std::cout << "GetValue_DFDX = " << element_i->GetValue(DFDX) << std::endl;
		}

		clock_t end = clock();
		std::cout << "  Objective Function sensitivities computed  [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		KRATOS_CATCH("");
	}


	/// Computes volume fraction sensitivities from the adjoint solution
	void ComputeVolumeFractionSensitivities_elem()
	{
		KRATOS_TRY;

		double Out = 0.0;

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			// Checking Calculate() - ok.
			element_i->Calculate(VOLUME, Out, ConstProcessInfo);
		}

		clock_t end = clock();
		std::cout << "  Volume fraction sensitivities computed     [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		KRATOS_CATCH("");
	}

		/// Computes compliance sensitivities from the adjoint solution (the problem is self-adjoint)
	std::vector<std::vector<double> > ComputeComplianceSensitivities_gpt()
	{
		KRATOS_TRY;

		std::vector<std::vector<double> > Sens_list;
		std::vector<double> Out; // dummy variables

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			// Checking CalculateOnIntegrationPoints() --  ok
			element_i->CalculateOnIntegrationPoints(COMPLIANCE, Out, ConstProcessInfo); 
			Sens_list.push_back(Out);
			// std::cout << "Out : " << Out << std::endl; 
		}


		clock_t end = clock();
		std::cout << "  Objective Function sensitivities computed  [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		return Sens_list;
		KRATOS_CATCH("");
	}


	/// Computes volume fraction sensitivities from the adjoint solution
	std::vector<std::vector<double> > ComputeVolumeFractionSensitivities_gpt()
	{
		KRATOS_TRY;

		std::vector<std::vector<double> > Sens_list;
		std::vector<double> Out; // dummy variables

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			// Checking Calculate() - ok.
			element_i->CalculateOnIntegrationPoints(VOLUME, Out, ConstProcessInfo);
			Sens_list.push_back(Out);
			// std::cout << Out << std::endl; 
		}

		clock_t end = clock();
		std::cout << "  Volume fraction sensitivities computed     [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		return Sens_list;
		KRATOS_CATCH("");
	}

	std::vector< std::vector< array_1d<double, 3 > > > ComputeGaussPointCoordinates()
	{
		KRATOS_TRY;

		std::vector< std::vector< array_1d<double, 3 > > > Gpts_lists;
		std::vector< array_1d<double, 3 > > Out; // dummy variables

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			// Checking Calculate() - ok.
			element_i-> CalculateOnIntegrationPoints(INTEGRATION_COORDINATES,  Out, ConstProcessInfo);
			Gpts_lists.push_back(Out);
			// std::cout << Out << std::endl; 
		}

		clock_t end = clock();
		std::cout << "  Gauss point coordinates computed     [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		return Gpts_lists;
		KRATOS_CATCH("");

	}


///@}

private:

	///@name Member Variables
	///@{

	ModelPart& mr_structure_model_part;
	ModelPart* mpAdjointModelPart;
	typename BaseType::Pointer mpStrategy;
	int m_dimension;

	///@}
	///@name Private Operations
	///@{

	///@}
}; // class StructureAdjointSensitivityStrategy

///@} // Kratos classes
///@} // AdjointStructureApplication group
}

#endif	/* KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED */
