// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Hayoung Chung
//
// ==============================================================================

#if !defined(KRATOS_RESIDUAL_BASED_STATIC_ERSATZ_SCHEME)
#define KRATOS_RESIDUAL_BASED_STATIC_ERSATZ_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

//#include "solid_mechanics_application.h"
#include "../../analytic_sensitivity_application.h"

// Application includes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */

/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class ResidualBasedIncrementalUpdateStaticErsatzScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedIncrementalUpdateStaticErsatzScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ResidualBasedIncrementalUpdateStaticErsatzScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>()
    {
    }

    /** Destructor.
    */
    virtual ~ResidualBasedIncrementalUpdateStaticErsatzScheme() {}

    /*@} */
    /**@name Operators
    */
    /*@{ */

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- UPDATE LHS AND RHS ----------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// This function calculates a new Youngs Modulus based on the densities and multiplies it into the
    /// LHS and RHS contributions of the complete system
    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType &LHS,
        LocalSystemVectorType &RHS,
        Element::EquationIdVectorType &EquationId,
        ProcessInfo &CurrentProcessInfo)
    {
        KRATOS_TRY
            //Initializing the non linear iteration for the current element
            (rCurrentElement)
                ->InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        //this computes LHS, RHS
        (rCurrentElement)->CalculateLocalSystem(LHS, RHS, CurrentProcessInfo);

        //Determine the new Youngs Modulus based on the assigned new density (X_PHYS)
        double rho_min = (rCurrentElement)->GetValue(RHO_MIN);
        double rho_current = (rCurrentElement)->GetValue(X_PHYS);
		if (rho_current <= rho_min )
			rho_current = rho_min ;
		if (rho_current > 1.0 || rho_current < 0.0)
			KRATOS_THROW_ERROR(std::logic_error, "current element density is outside [0, 1] ", rho_current);
        double factor = rho_current; 

        // Factorize LHS and RHS according SIMP approach
        // Note that when this function is called, all the contributions from the force conditions are missing.
        // I.e. RHS = -K*u_init. Hence we can directly factorize LHS and RHS to obtained the modified stiffnesses
        LHS *= factor;
        RHS *= factor;

        //Continuation of the basic operations
        (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /*@} */
    /**@name Operations */
    /*@{ */

    /*@} */
    /**@name Access */
    /*@{ */

    /*@} */
    /**@name Inquiry */
    /*@{ */

    /*@} */
    /**@name Friends */
    /*@{ */

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /*@} */
    /**@name Private  Access */
    /*@{ */

    /*@} */
    /**@name Private Inquiry */
    /*@{ */

    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_STATIC_ERSATZ_SCHEME  defined */
