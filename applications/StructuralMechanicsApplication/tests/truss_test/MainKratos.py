from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis as sld

import time

if __name__ == "__main__":
    with open("nonlinear_3D2NTruss_plastic_snapthrough_test_parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = sld.StructuralMechanicsAnalysis(model,parameters)
    simulation.Run()