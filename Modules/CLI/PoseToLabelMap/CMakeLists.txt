
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
  TARGET_LIBRARIES ${ITK_LIBRARIES} vtkIO vtkGraphics ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
