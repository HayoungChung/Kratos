set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosAnalyticSensitivityApplication ****")

include(pybind11Tools)


include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../ConstitutiveModelsApplication/ )
include_directories( /home/hayoung/anaconda3/envs/kratos/include/python3.7m/)


#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication/ )

## generate variables with the sources
file(GLOB_RECURSE KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/analytic_sensitivity_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_ersatz_element.cpp
    # ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    # file(GLOB_RECURSE KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

##  python interface sources
file(GLOB_RECURSE KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_PYTHON_INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/analytic_sensitivity_python_application.cpp)

add_library(KratosAnalyticSensitivityCore SHARED ${KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_CORE} )# ${KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosAnalyticSensitivityCore PUBLIC KratosCore KratosSolidMechanicsCore)
set_target_properties(KratosAnalyticSensitivityCore PROPERTIES COMPILE_DEFINITIONS "ANALYTIC_SENSITIVITY_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosAnalyticSensitivityApplication MODULE THIN_LTO ${KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosAnalyticSensitivityApplication PUBLIC KratosAnalyticSensitivityCore KratosSolidMechanicsCore)
set_target_properties(KratosAnalyticSensitivityApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosAnalyticSensitivityApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosAnalyticSensitivityApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    cotire(KratosAnalyticSensitivityCore)
    cotire(KratosAnalyticSensitivityApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosAnalyticSensitivityCore DESTINATION libs )
install(TARGETS KratosAnalyticSensitivityApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/AnalyticSensitivityApplication.py" DESTINATION "KratosMultiphysics/AnalyticSensitivityApplication" RENAME "__init__.py")


# ###
# set( KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_SOURCES
# 	${CMAKE_CURRENT_SOURCE_DIR}/analytic_sensitivity_application.cpp

# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/analytic_sensitivity_python_application.cpp
# 	# ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
# 	# ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_ersatz_element.cpp
# )



# ###############################################################
# ## define library Kratos which defines the basic python interface
# add_library(KratosAnalyticSensitivityApplication SHARED ${KRATOS_ANALYTIC_SENSITIVITY_APPLICATION_SOURCES})
# target_link_libraries(KratosAnalyticSensitivityApplication KratosCore KratosSolidMechanicsApplication)
# set_target_properties(KratosAnalyticSensitivityApplication PROPERTIES PREFIX "")
# install(TARGETS KratosAnalyticSensitivityApplication DESTINATION libs )


# # get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# # message("TestApplication subdir inc_dirs = ${inc_dirs}")

# # changing the .dll suffix to .pyd
# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# 	set_target_properties(KratosAnalyticSensitivityApplication PROPERTIES SUFFIX .pyd)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# if(${INSTALL_PYTHON_FILES} MATCHES ON)
#   get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#   install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
# endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# # get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# # message("TestApplication subdir inc_dirs = ${inc_dirs}")

# # Add to the KratosMultiphisics Python module
# install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/AnalyticSensitivityApplication.py" DESTINATION KratosMultiphysics )
