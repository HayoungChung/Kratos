# ==============================================================================
#  TopologyOptimizationApplication
#
#  License:         BSD License
#                   license: AnalyticSensitivityApplication/license.txt
#
#  Main authors:    Hayoung Chung
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
from KratosAnalyticSensitivityApplication import *
application = KratosAnalyticSensitivityApplication()
application_name = "KratosAnalyticSensitivityApplication"
application_folder = "AnalyticSensitivityApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)