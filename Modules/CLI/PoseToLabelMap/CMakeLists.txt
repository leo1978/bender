
#-----------------------------------------------------------------------------
set(MODULE_NAME PoseToLabelMap) # Do not use 'project()'

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  ADDITIONAL_SRCS vtkPolyDataPointIDSampler.cxx
                  vtkPolyDataPointIDSampler.h
                  vtkPolyDataPointSampler.cxx
                  vtkPolyDataPointSampler.h
  INCLUDE_DIRECTORIES ${vtkSlicerArmaturesModuleLogic_SOURCE_DIR}
                      ${vtkSlicerArmaturesModuleLogic_BINARY_DIR}
                      ${vtkSlicerModelsModuleLogic_SOURCE_DIR}
                      ${vtkSlicerModelsModuleLogic_BINARY_DIR}
                      ${SlicerBaseLogic_SOURCE_DIR}
                      ${SlicerBaseLogic_BINARY_DIR}
                      ${MRMLCore_SOURCE_DIR}
                      ${MRMLCore_BINARY_DIR}
                      ${MRMLLogic_SOURCE_DIR}
                      ${MRMLLogic_BINARY_DIR}
  TARGET_LIBRARIES ${ITK_LIBRARIES} vtkIO vtkGraphics
                   vtkSlicerArmaturesModuleLogic
                   vtkSlicerModelsModuleLogic
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
