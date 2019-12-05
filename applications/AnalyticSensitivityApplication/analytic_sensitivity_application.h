// ==============================================================================
//  KratosAnalyticSensitivityApplication
//
//  License:         BSD License
//                   license: AnalyticSensitivityApplication/license.txt
//
//  Main authors:   Hayoung Chung,  
//                  Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#if !defined(KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_H_INCLUDED )
#define  KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Core applications
// #include "analytic_sensitivity_application_variables.h"
// #include "analytic_sensitivity_application.h"

#include "./custom_elements/small_displacement_ersatz_element.hpp"

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

//elements


namespace Kratos
{

	///@name Kratos Globals
	///@{

    // Variables definition with Python connection
	KRATOS_DEFINE_VARIABLE( double, RHO_MIN )
    KRATOS_DEFINE_VARIABLE( double, X_PHYS ) // material density
    KRATOS_DEFINE_VARIABLE( double, X_PHYS_OLD ) 
    KRATOS_DEFINE_VARIABLE( double, DFDX ) // sensitivity in an element
    KRATOS_DEFINE_VARIABLE( Vector, DFDX_intg ) // sensitivity in an element
    KRATOS_DEFINE_VARIABLE( double, SOLID_VOID )
    KRATOS_DEFINE_VARIABLE( double, LOCAL_STRAIN_ENERGY )
    KRATOS_DEFINE_VARIABLE( double, COMPLIANCE )
	KRATOS_DEFINE_VARIABLE( Vector, COMPLIANCE_intg )
	KRATOS_DEFINE_VARIABLE( double, VOLUME )
	KRATOS_DEFINE_VARIABLE( double, P_NORM_STRESS )
	KRATOS_DEFINE_VARIABLE( double, UF )


	///@}
	///@name Type Definitions
	///@{

	///@}
	///@name  Enum's
	///@{

	///@}
	///@name  Functions
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KratosAnalyticSensitivityApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{


		/// Pointer definition of KratosAnalyticSensitivityApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosAnalyticSensitivityApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosAnalyticSensitivityApplication();

		/// Destructor.
		~KratosAnalyticSensitivityApplication() override {}


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



		///@}
		///@name Access
		///@{


		///@}
		///@name Inquiry
		///@{


		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosAnalyticSensitivityApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
      	KRATOS_WATCH("in my application");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
      }


		///@}
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables
		///@{


		///@}
		///@name Protected member Variables
		///@{


		///@}
		///@name Protected Operators
		///@{


		///@}
		///@name Protected Operations
		///@{


		///@}
		///@name Protected  Access
		///@{


		///@}
		///@name Protected Inquiry
		///@{


		///@}
		///@name Protected LifeCycle
		///@{


		///@}

	private:
		///@name Static Member Variables
		///@{

		///@}
		///@name Member Variables
		///@{

        //small_displacement
      	const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D3N; // dummy element for surface representation
        const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D4N;
        const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D8N;
		const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement2D4N;

//        Extra elements to be added in the future
//        const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D6N;
//        const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D10N;
//        const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D15N;
//        const SmallDisplacementErsatzElement mSmallDisplacementErsatzElement3D20N;
//        const SmallDisplacementErsatzElement mSmallDis3D27N;


		///@}
		///@name Private Operators
		///@{


		///@}
		///@name Private Operations
		///@{


		///@}
		///@name Private  Access
		///@{


		///@}
		///@name Private Inquiry
		///@{


		///@}
		///@name Un accessible methods
		///@{

		/// Assignment operator.
		KratosAnalyticSensitivityApplication& operator=(KratosAnalyticSensitivityApplication const& rOther);

		/// Copy constructor.
		KratosAnalyticSensitivityApplication(KratosAnalyticSensitivityApplication const& rOther);


		///@}

	}; // Class KratosAnalyticSensitivityApplication

	///@}


	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{

	///@}


}  // namespace Kratos.

#endif // KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_H_INCLUDED  defined


