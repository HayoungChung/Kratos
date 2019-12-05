from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructure
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers
import KratosMultiphysics.AnalyticSensitivityApplication as KratosSensitivity# without this, Ersatz element cannot be used
import structural_mechanics_static_Ersatz_solver
import numpy as np
import itertools
import pickle

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

class SensitivityField:
    def __init__(self, Gpts_coord, nCons = 1):
        self.nELEM = len(Gpts_coord)
        self.nOrder = len(Gpts_coord[0])
        self.nGpts = self.nELEM * self.nOrder
        self.Sf = []
        self.Sg = [None]*nCons
        self.Coords = list(itertools.chain.from_iterable(Gpts_coord))
    def setSensObj(self, Sf):
        if len(Sf)*len(Sf[0]) != self.nGpts:
            print("nGpts error")
            raise
        self.Sf = Sf
    def setSensCons(self, Sg, index = 0):
        if len(Sg)*len(Sg[0]) != self.nGpts:
            print("nGpts error")
            raise
        self.Sg[index] = Sg
    
    

class ErsatzAnalysis(StructuralMechanicsAnalysis):
    """ 
    This class is the main script that runs ErsatzAnalysis
    """
    def __init__(self, model, project_parameters):
        super(ErsatzAnalysis, self).__init__(model, project_parameters) 
        self.model = model
        self.project_parameters

    def _CreateSolver(self):
        solver = structural_mechanics_static_Ersatz_solver.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetSolutionStrategy(self):
        linear_solver = self._GetSolver().get_linear_solver(); # linear_solver_factory.ConstructSolver(ProjectParameters["solver_settings"]["linear_solver_settings"])
        self.ErsatzStrategy = KratosSensitivity.StructureAdjointSensitivityStrategy(self.model.GetModelPart('Structure'), linear_solver, 3)
        return self.ErsatzStrategy

    # def _CreateAdjointSolver(self):
    #     solver =

if __name__ == "__main__":
    from sys import argv

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ErsatzAnalysis(model, parameters)
    simulation.Initialize()

    model_part = model.GetModelPart('Structure')
    elem_vec = model_part.GetElements()
    for element_i in elem_vec:
            element_i.SetValue(KratosMultiphysics.AnalyticSensitivityApplication.RHO_MIN, 1e-4)
            element_i.SetValue(KratosMultiphysics.AnalyticSensitivityApplication.X_PHYS, 0.01)
            element_i.SetValue(KratosMultiphysics.AnalyticSensitivityApplication.X_PHYS_OLD, 0.7)

    simulation.RunSolutionLoop()

    SS = simulation._GetSolver().get_solution_scheme()
    TT = simulation._GetSolutionStrategy()

    TT.ComputeComplianceSensitivities_elem()
    TT.ComputeVolumeFractionSensitivities_elem()

    Sf = TT.ComputeComplianceSensitivities_gpt()
    Sg = TT.ComputeVolumeFractionSensitivities_gpt()
    Gpts_list = TT.ComputeGaussPointCoordinates()

    AA = SensitivityField(Gpts_list)
    AA.setSensCons(Sg)
    AA.setSensObj(Sf)


    for ii in range(AA.nELEM):
        print(AA.Coords[ii])
        print(AA.Sf[ii])
        print(AA.Sg[0][ii])
    
    # print(Gpts_list)

    # np.savetxt("Sf.txt", Sf)
    # np.savetxt("Sg.txt", Sg)
    # np.savetxt("GG.txt", Gpts_list)



    # # print(simulation._GetSolver()._solution_strategy)
    # # help(simulation._GetSolver())
    # elem_vec = model_part.GetElements()

    ## PRINT!!
    from KratosMultiphysics.vtk_output_process import VtkOutputProcess
    print(parameters["output_processes"]["vtk_output"][0]["Parameters"]["model_part_name"].GetString())
    vtkprocess = VtkOutputProcess(model, parameters["output_processes"]["vtk_output"][0]["Parameters"])
    vtkprocess.ExecuteInitialize()
    vtkprocess.ExecuteBeforeSolutionLoop()
    vtkprocess.ExecuteInitializeSolutionStep()
    vtkprocess.PrintOutput()
    vtkprocess.ExecuteFinalizeSolutionStep()
    vtkprocess.ExecuteFinalize()


    vtkprocess2 = VtkOutputProcess(model, KratosMultiphysics.Parameters("""{
                                                "model_part_name"                    : "Structure",
                                                "gauss_point_variables_in_elements": ["COMPLIANCE"],
                                                "folder_name"                        : "vtk_output2"
                                            }
                                            """))
    vtkprocess2.ExecuteInitialize()
    vtkprocess2.ExecuteBeforeSolutionLoop()
    vtkprocess2.ExecuteInitializeSolutionStep()
    vtkprocess2.PrintOutput()
    vtkprocess2.ExecuteFinalizeSolutionStep()
    vtkprocess2.ExecuteFinalize()


    # for ee in elem_vec:
    #     # print(ee.GetIntegrationPoints()) # works ok
    #     # print(ee.GetValuesOnIntegrationPoints(KratosMultiphysics.AnalyticSensitivityApplication.COMPLIANCE, model_part.ProcessInfo))
    #     # TODO: BELOW DO NOT WORK.
    #     print(ee.CalculateOnIntegrationPoints(KratosStructure.VON_MISES_STRESS, model_part.ProcessInfo)) # zero as well
    #     print(ee.CalculateOnIntegrationPoints(KratosSensitivity.COMPLIANCE, model_part.ProcessInfo)) # 
    #     print(ee.CalculateOnIntegrationPoints(KratosSensitivity.COMPLIANCE_intg, model_part.ProcessInfo)[0])
    #     print(ee.CalculateOnIntegrationPoints(KratosSensitivity.COMPLIANCE_intg, model_part.ProcessInfo)[1])
    #     print(ee.CalculateOnIntegrationPoints(KratosSensitivity.COMPLIANCE_intg, model_part.ProcessInfo)[2])
    #     print(ee.CalculateOnIntegrationPoints(KratosSensitivity.COMPLIANCE_intg, model_part.ProcessInfo)[3])
    #     print(ee.CalculateOnIntegrationPoints(KratosSensitivity.VOLUME, model_part.ProcessInfo)[0])
    #     AssertionError("FIXME: these errors possibly originated from 2D-3D element compatibility")


    #     # print(ee.GetValuesOnIntegrationPoints(KratosSensitivity.COMPLIANCE_intg, model_part.ProcessInfo))

    #     # BELOW WORKS OK.
    #     print(ee.GetValue(KratosSensitivity.COMPLIANCE))
    #     print(ee.GetValue(KratosSensitivity.VOLUME))
    #     # print(ee.GetValuesOnIntegrationPoints(KratosStructure.VON_MISES_STRESS, model_part.ProcessInfo))
    #     # print(ee.GetValue(KratosMultiphysics.AnalyticSensitivityApplication.COMPLIANCE))
    #     # print(ee.GetValue(KratosMultiphysics.AnalyticSensitivityApplication.X_PHYS))
    #     pass

    # node_vec = model_part.GetNodes()
    # for nn in node_vec:
    #     pass
    #     # print(nn.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)) # this works. 
    #     # print(nn.GetValue(KratosMultiphysics.DISPLACEMENT)) # this does not.

    simulation.Finalize()
    
