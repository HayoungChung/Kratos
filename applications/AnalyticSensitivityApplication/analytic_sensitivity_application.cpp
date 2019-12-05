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

// System includes

// External includes

// Project includes
#include <vector>

#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/variables.h"

#include "./analytic_sensitivity_application.h"

// Geometries that must be added when more elements are added into the application (SOLID MECHANICS APPLICATION)
//#include "geometries/tetrahedra_3d_10.h"
//#include "geometries/hexahedra_3d_20.h"
//#include "geometries/hexahedra_3d_27.h"
//#include "geometries/prism_3d_6.h"
//#include "geometries/prism_3d_15.h"

namespace Kratos
{
// Variables definition with Python connection
KRATOS_CREATE_VARIABLE(double, RHO_MIN)
KRATOS_CREATE_VARIABLE(double, X_PHYS) // material density
KRATOS_CREATE_VARIABLE(double, X_PHYS_OLD)
KRATOS_CREATE_VARIABLE(double, DFDX)                   // sensitivity in an element
KRATOS_CREATE_VARIABLE(Vector, DFDX_intg) // sensitivity in an element
KRATOS_CREATE_VARIABLE(double, SOLID_VOID)
KRATOS_CREATE_VARIABLE(double, LOCAL_STRAIN_ENERGY)
KRATOS_CREATE_VARIABLE(double, COMPLIANCE)
KRATOS_CREATE_VARIABLE(Vector, COMPLIANCE_intg)
KRATOS_CREATE_VARIABLE(double, VOLUME)
KRATOS_CREATE_VARIABLE(double, P_NORM_STRESS)
KRATOS_CREATE_VARIABLE(double, UF)

KratosAnalyticSensitivityApplication::KratosAnalyticSensitivityApplication() : KratosApplication("AnalyticSensitivityApplication"),
                                                                               mSmallDisplacementErsatzElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))), // dummy element for surface representation
                                                                               mSmallDisplacementErsatzElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
                                                                               mSmallDisplacementErsatzElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
                                                                               mSmallDisplacementErsatzElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))

//        Extra elements that can be added in the future
//        mSmallDisplacementErsatzElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
//        mSmallDisplacementErsatzElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
//        mSmallDisplacementErsatzElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
//        mSmallDisplacementErsatzElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
//        mSmallDisplacementErsatzElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )

{
}

void KratosAnalyticSensitivityApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::cout << "Initializing KratosAnalyticSensitivityApplication...    " << std::endl;

    /* ==================================================================================================== 
        Register small displacement elements
    ===================================================================================================== */
    KRATOS_REGISTER_ELEMENT("SmallDisplacementErsatzElement3D3N", mSmallDisplacementErsatzElement3D3N) // dummy element for surface representation
    KRATOS_REGISTER_ELEMENT("SmallDisplacementErsatzElement3D4N", mSmallDisplacementErsatzElement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementErsatzElement3D8N", mSmallDisplacementErsatzElement3D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementErsatzElement2D4N", mSmallDisplacementErsatzElement2D4N)

    //        Extra elements that can be added in the future
    //        KRATOS_REGISTER_ELEMENT( "SmallDisplacementErsatzElement3D6N", mSmallDisplacementErsatzElement3D6N )
    //        KRATOS_REGISTER_ELEMENT( "SmallDisplacementErsatzElement3D10N", mSmallDisplacementErsatzElement3D10N )
    //        KRATOS_REGISTER_ELEMENT( "SmallDisplacementErsatzElement3D15N", mSmallDisplacementErsatzElement3D15N )
    //        KRATOS_REGISTER_ELEMENT( "SmallDisplacementErsatzElement3D20N", mSmallDisplacementErsatzElement3D20N )
    //        KRATOS_REGISTER_ELEMENT( "SmallDisplacementErsatzElement3D27N", mSmallDisplacementErsatzElement3D27N )

    /* ==================================================================================================== 
        Register Variables with Python connection
    ===================================================================================================== */
    KRATOS_REGISTER_VARIABLE(RHO_MIN)
    KRATOS_REGISTER_VARIABLE(X_PHYS) // material density
    KRATOS_REGISTER_VARIABLE(X_PHYS_OLD)
    KRATOS_REGISTER_VARIABLE(DFDX) // sensitivity in an element
    KRATOS_REGISTER_VARIABLE(DFDX_intg) // sensitivity in an element
    KRATOS_REGISTER_VARIABLE(SOLID_VOID)
    KRATOS_REGISTER_VARIABLE(LOCAL_STRAIN_ENERGY)
    KRATOS_REGISTER_VARIABLE(COMPLIANCE)
    KRATOS_REGISTER_VARIABLE(COMPLIANCE_intg)
    KRATOS_REGISTER_VARIABLE(VOLUME)
    KRATOS_REGISTER_VARIABLE(P_NORM_STRESS)
    KRATOS_REGISTER_VARIABLE(UF)
}

} // namespace Kratos.
