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

// Application includes
#include "small_displacement_ersatz_element.hpp"

#include "../analytic_sensitivity_application.h"

namespace Kratos
{
// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================
SmallDisplacementErsatzElement ::SmallDisplacementErsatzElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: SmallDisplacementElement(NewId, pGeometry)
{
}

SmallDisplacementErsatzElement ::SmallDisplacementErsatzElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: SmallDisplacementElement(NewId, pGeometry, pProperties)
{
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

// =============================================================================================================================================
// COPY CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementErsatzElement::SmallDisplacementErsatzElement(SmallDisplacementErsatzElement const &rOther)
	: SmallDisplacementElement(rOther)
{
}

// =============================================================================================================================================
// OPERATIONS
// =============================================================================================================================================

Element::Pointer SmallDisplacementErsatzElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	return Element::Pointer(new SmallDisplacementErsatzElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

// =============================================================================================================================================
// CLONE
// =============================================================================================================================================

Element::Pointer SmallDisplacementErsatzElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

	SmallDisplacementErsatzElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

	//-----------//

	NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

	if (NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size())
	{
		NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

		if (NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber())
			KRATOS_THROW_ERROR(std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size());
	}

	return Element::Pointer(new SmallDisplacementErsatzElement(NewElement));
}

// =============================================================================================================================================
// DESTRUCTOR
// =============================================================================================================================================

SmallDisplacementErsatzElement::~SmallDisplacementErsatzElement()
{
}

// =============================================================================================================================================
// STARTING / ENDING METHODS
// =============================================================================================================================================

//************************************************************************************
//************************************************************************************
void SmallDisplacementErsatzElement::GetValueOnIntegrationPoints(const Variable<double> &rVariable,
																 std::vector<double> &rValues,
																 const ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	// Additional part for post-processing of the topology optimized model part
	if (rVariable == X_PHYS)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	else if (rVariable == COMPLIANCE)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	else if (rVariable == COMPLIANCE_intg)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	else if (rVariable == VOLUME)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);

	// From original SmallDisplacementElement
	else if (rVariable == VON_MISES_STRESS)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	else
	{

		const unsigned int &integration_points_number = GetGeometry()
															.IntegrationPointsNumber(mThisIntegrationMethod);

		if (rValues.size() != integration_points_number)
			rValues.resize(integration_points_number);

		for (unsigned int ii = 0; ii < integration_points_number; ii++)
		{
			rValues[ii] = 0.0;
			rValues[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable,
															   rValues[ii]);
		}
	}

	KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementErsatzElement::Calculate(const Variable<double> &rVariable, double &rOutput, const ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	if (rVariable == LOCAL_STRAIN_ENERGY || rVariable == COMPLIANCE)
	{
		// Get values
		double rho_min = this->GetValue(RHO_MIN );
		if (rho_min >= 0.1) 
			KRATOS_THROW_ERROR(std::logic_error, "Minimum value of area is too large ", rho_min );

		// for compliance minimization
		MatrixType Ke0 = Matrix();
		this->CalculateLeftHandSide(Ke0, const_cast<ProcessInfo &>(rCurrentProcessInfo));
		double rho_current = this->GetValue(X_PHYS);

		if (rho_current <= rho_min )
			rho_current = rho_min ;
		if (rho_current > 1.0 || rho_current < 0.0)
			KRATOS_THROW_ERROR(std::logic_error, "current element density is outside [0, 1] ", rho_current);
		MatrixType Ke = Ke0 * rho_current;

		// Loop through nodes of elements and create elemental displacement vector "ue"
		Element::GeometryType &rGeom = this->GetGeometry();
		unsigned int NumNodes = rGeom.PointsNumber(); //NumNodes=8

		// Resize "ue" according to element type
		Vector ue;
		ue.resize(NumNodes * 2);

		// Set the displacement obtained from the FE-Analysis
		for (unsigned int node_i = 0; node_i < NumNodes; node_i++)
		{
			array_1d<double, 3> &CurrentDisplacement = rGeom[node_i].FastGetSolutionStepValue(DISPLACEMENT);
			ue[2 * node_i + 0] = CurrentDisplacement[0];
			ue[2 * node_i + 1] = CurrentDisplacement[1];
			// ue[2 * node_i + 2] = CurrentDisplacement[2];
		}

		// Calculate trans(ue)*Ke0*ue
		Vector intermediateVector;
		intermediateVector.resize(NumNodes * 3);
		intermediateVector = prod(trans(ue), Ke);
		double ue_Ke_ue = inner_prod(intermediateVector, ue);

		// Calculation of the compliance sensitivities DCDX
		double dfdx = (-1.0) * ue_Ke_ue;
		double local_strain_energy = 0.5 * ue_Ke_ue;
		// Calculation of the local strain energy (requires Ke)
		if (rVariable == COMPLIANCE)
		{
			rOutput = dfdx;
			this->SetValue(DFDX, dfdx);
			this->SetValue(COMPLIANCE, dfdx);
		}
		else if (rVariable == LOCAL_STRAIN_ENERGY)
		{
			rOutput = local_strain_energy;
			this->SetValue(LOCAL_STRAIN_ENERGY, local_strain_energy);
			this->SetValue(COMPLIANCE, local_strain_energy);
		}
	}
	else if (rVariable == VOLUME)
	{
		// Calculation of the volume sensitivities DVDX
		this->SetValue(DFDX, 1.0);
		double x_phys = this->GetValue(X_PHYS);
		this->SetValue(VOLUME, x_phys);
		rOutput = x_phys;
	}

	KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementErsatzElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    // if (mConstitutiveLawVector[0]->Has( rVariable)) {
    //     GetValueOnConstitutiveLaw(rVariable, rOutput);
    // } else {
        if (rVariable == INTEGRATION_COORDINATES) {
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            // KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                Point global_point;
                GetGeometry().GlobalCoordinates(global_point, integration_points[point_number]);

                rOutput[point_number] = global_point.Coordinates();
            }
        } //else {
            // CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        // }
    // }
}

void SmallDisplacementErsatzElement::CalculateOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rOutput,
																  const ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	// From original SmallDisplacementElement
	const unsigned int &integration_points_number = GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

	if (rOutput.size() != integration_points_number)
		rOutput.resize(integration_points_number, false);

	if (rVariable == VON_MISES_STRESS)
	{
		//create and initialize element variables:
		ElementDataType Variables;
		InitializeElementData(Variables, rCurrentProcessInfo);

		//create constitutive law parameters:
		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(),
										   rCurrentProcessInfo);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();

		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

		for (unsigned int PointNumber = 0;
			 PointNumber < mConstitutiveLawVector.size(); PointNumber++)
		{
			//compute element kinematics B, F, DN_DX ...
			CalculateKinematics(Variables, PointNumber);

			//set general variables to constitutivelaw parameters
			SetElementData(Variables, Values, PointNumber);

			//call the constitutive law to update material variables
			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

			// ComparisonUtilities EquivalentStress;
			// rOutput[PointNumber] = EquivalentStress.CalculateVonMises(
			// 	Variables.StressVector);
			rOutput[PointNumber] =  ElementUtilities::CalculateVonMises(Variables.StressVector);
		}
		// this->SetValueOnIntegrationPoints(VON_MISES_STRESS, rOutput, rCurrentProcessInfo); // setValue is still virtual: not defined
	}

	// Additional part for post-processing of the topology optimized model part
	else if (rVariable == X_PHYS)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
			rOutput[PointNumber] = this->GetValue(X_PHYS);
	}
	// Additonal part for sensitivities evaluated at the gausspoints (Hayoung)
	else if (rVariable == VOLUME)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
		{
			rOutput[PointNumber] = 1.0;
		}
		// this->SetValueOnIntegrationPoints(VOLUME, rOutput, rCurrentProcessInfo); // setValue is still virtual: not defined

	}
	else if (rVariable == COMPLIANCE)
	{
		//create and initialize element variables:
		ElementDataType Variables;
		this->InitializeElementData(Variables, rCurrentProcessInfo);

		//create constitutive law parameters (Initialize pointers to NULL)
		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(),
										   rCurrentProcessInfo);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();


		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
		ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

		for (unsigned int PointNumber = 0;
			 PointNumber < mConstitutiveLawVector.size(); PointNumber++)
		{
			//compute element kinematics B, F, DN_DX ...
			this->CalculateKinematics(Variables, PointNumber);

			//set general variables to constitutivelaw parameters
			this->SetElementData(Variables, Values, PointNumber);

			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values); // based on detF (which is one.)

			// it works.
			// std::cout << "Variables.StressVector: " << Variables.StressVector << std::endl;
			// std::cout << "Variables.StrainVector: " << Variables.StrainVector << std::endl; 

			double localE = inner_prod(Variables.StressVector, Variables.StrainVector);
			rOutput[PointNumber] = -localE;			
		}
		// this->SetValueOnIntegrationPoints(COMPLIANCE, rOutput, rCurrentProcessInfo); // setValue is still virtual: not defined

	}
	else
	{
		for (unsigned int ii = 0; ii < integration_points_number; ii++)
			rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rOutput[ii]);
	}
	KRATOS_CATCH("")
}
// =============================================================================================================================================
// =============================================================================================================================================

void SmallDisplacementErsatzElement::save(Serializer &rSerializer) const
{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementElement)
}

void SmallDisplacementErsatzElement::load(Serializer &rSerializer)
{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementElement)
}

} // Namespace Kratos