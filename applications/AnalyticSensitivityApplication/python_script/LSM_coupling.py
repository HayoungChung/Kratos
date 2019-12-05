from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructure
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers
import KratosMultiphysics.AnalyticSensitivityApplication as KratosSensitivity# without this, Ersatz element cannot be used
import structural_mechanics_static_Ersatz_solver
import numpy as np
import itertools
import pickle
import pyslsm as lsm # 2d level-set method library
import matplotlib.pyplot as plt
from scipy import interpolate 

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

    level_set = lsm.LevelSet(160, 80)
    mesh = level_set.mesh
    boundary = lsm.Boundary()
    vtk_io = lsm.InputOutput()
    boundary.discretise(level_set)
    # print(boundary.length)
    # boundary.ComputeNormalVectors(level_set)
    vtk_io.saveLevelSetVTK("aa_1.vtk", level_set)
    level_set.computeAreaFractions(boundary)
    # vtk_io.saveAreaFractionsVTK("bb_1.vtk", level_set.mesh)
    # for ii in level_set.mesh.elements:
    #     print(ii.area)
    

    model = KratosMultiphysics.Model()
    simulation = ErsatzAnalysis(model, parameters)
    simulation.Initialize()

    model_part = model.GetModelPart('Structure')
    elem_vec = model_part.GetElements()
    cnt = 0 
    for element_i in elem_vec:
        element_i.SetValue(KratosMultiphysics.AnalyticSensitivityApplication.RHO_MIN, 1e-4)
        element_i.SetValue(KratosMultiphysics.AnalyticSensitivityApplication.X_PHYS, level_set.mesh.elements[cnt].area)
        element_i.SetValue(KratosMultiphysics.AnalyticSensitivityApplication.X_PHYS_OLD, 0.7)
        cnt += 1

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

    # interpolation ===========
    xx = []
    yy = []
    sfs = []
    sgs = []

    for ii in range(AA.nGpts):
        xx.append(AA.Coords[ii][0])
        yy.append(AA.Coords[ii][1])
    for ii in range(AA.nELEM):
        for jj in range(AA.nOrder):
            # xx.append(AA.Coords[ii][jj])
            # yy.append(AA.Coords[ii][jj])
            sfs.append(AA.Sf[ii][jj])
            sgs.append(AA.Sg[0][ii][jj])

    F_field = interpolate.SmoothBivariateSpline(xx,yy,sfs)
    G_field = interpolate.SmoothBivariateSpline(xx,yy,sgs)

    boundary_points = boundary.points
    bnd_xx = []
    bnd_yy = []
    for ii in boundary_points:
        bnd_xx.append(ii.coord.x)
        bnd_yy.append(ii.coord.y)
    
    plt.figure()
    plt.scatter(xx,yy,10,sfs,marker='o')
    plt.scatter(bnd_xx,bnd_yy,10,F_field.ev(bnd_xx,bnd_yy),marker='x')

    # plt.figure()
    # plt.scatter(xx,yy,10,sgs)
    plt.show()

    
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


    simulation.Finalize()
    
