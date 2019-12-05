from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
from KratosMultiphysics.AnalyticSensitivityApplication import ResidualBasedIncrementalUpdateStaticErsatzScheme, StructureAdjointSensitivityStrategy

def CreateSolver(model, custom_settings):
    return StaticMechanicalErsatzSolver(model, custom_settings)

class StaticMechanicalErsatzSolver(MechanicalSolver):
    """The structural mechanics static solver using Ersatz material

    This class creates the mechanical solvers for static analysis.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(StaticMechanicalErsatzSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticMechanicalErsatzSolver]:: ", "Construction finished")

    def _create_solution_scheme(self):
        return ResidualBasedIncrementalUpdateStaticErsatzScheme()

    # def _create_linear_strategy(self):
    #     computing_model_part = self.GetComputingModelPart()
    #     mechanical_scheme = self.get_solution_scheme()
    #     linear_solver = self.get_linear_solver()
    #     builder_and_solver = self.get_builder_and_solver()
    #     return StructureAdjointSensitivityStrategy(computing_model_part, linear_solver, 3)
    #                                                         #   builder_and_solver,
    #                                                         #   self.settings["compute_reactions"].GetBool(),
    #                                                         #   self.settings["reform_dofs_at_each_step"].GetBool(),
    #                                                         #   False,
    #                                                         #   self.settings["move_mesh_flag"].GetBool())