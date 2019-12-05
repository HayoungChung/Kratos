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

#if !defined(KRATOS_SMALL_DISPLACEMENT_ERSTAZ_ELEMENT_H_INCLUDED)
#define KRATOS_SMALL_DISPLACEMENT_ERSTAZ_ELEMENT_H_INCLUDED

// project includes
// #include "includes/define.h"
#include "/home/hayoung/packages/Kratos/kratos/includes/define.h"
#include "/home/hayoung/packages/Kratos/applications/SolidMechanicsApplication/solid_mechanics_application.h"
#include "/home/hayoung/packages/Kratos/applications/SolidMechanicsApplication/custom_elements/solid_elements/small_displacement_element.hpp"
// #include "/home/hayoung/packages/Kratos/applications/StructuralMechanicsApplication/structural_mechanics_application.h"
// #include "/home/hayoung/packages/Kratos/applications/StructuralMechanicsApplication/custom_elements/small_displacement.h"

namespace Kratos
{
///@name Kratos Globals
///@{
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

/// Sensitivity analysis module that uses small displacement elements 
// this substitutes adjoint sensitivity and topology optimization modules

/**
 * implements small-displacement lagrangian element for structural analysis
 * this should work for arbitrary geometries in 2d and 3d, but currently only Q4 mesh is implemented 
 * #TODO: check if this works for unstructured mesh
 */
class KRATOS_API(ANALYTIC_SENSITIVITY_APPLICATION) SmallDisplacementErsatzElement : public SmallDisplacementElement
{
    public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    ///Type for element variables
    typedef SmallDisplacementElement::ElementDataType ElementDataType;

    /// Counted pointer of SmallDisplacementErsatzElement

    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementErsatzElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    SmallDisplacementErsatzElement() : SmallDisplacementElement()
    {};
    
    SmallDisplacementErsatzElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SmallDisplacementErsatzElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    SmallDisplacementErsatzElement(SmallDisplacementErsatzElement const& rOther);

    /// Destructor.
    ~SmallDisplacementErsatzElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator (overload)
    SmallDisplacementErsatzElement& operator=(SmallDisplacementErsatzElement const& rOther);
    
    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     * TODO: not implemented yet
     */
    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

    // =============================================================================================================================================
    // STARTING / ENDING METHODS
    // =============================================================================================================================================

    /**
     * Function that gets the value on the Integration Point (For printing purposes in the output)
     * @param rVariable: variable type flag
     * @param rValues: return values
     * @param rCurrentProcessInfo : Process
     */ 
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Function to calculate the sensitivities and the objective function
     * @param rVariable: variable type flag
     * @param rValues: return values
     * @param rCurrentProcessInfo : Process
     */ 
    void Calculate(const Variable<double> &rVariable, double &rOutput, const ProcessInfo &rCurrentProcessInfo);

    
    /**
     * Function that overwrites the CalculateOnIntegrationPoints, to insert the X_PHYS into all Gauss Points of the given element
     * That allows printing X_PHYS as elemental value in GiD
     * @param rVariable: variable type flag
     * @param rValues: return values
     * @param rCurrentProcessInfo : Process
     */ 
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
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

    // SmallDisplacementErsatzElement() : SmallDisplacementElement()
    // {};

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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    virtual void save(Serializer& rSerializer) const override;
    virtual void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class SmallDisplacementSIMPElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif