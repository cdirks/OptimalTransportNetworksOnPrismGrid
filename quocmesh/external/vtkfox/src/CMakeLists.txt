SET(VTKFOX_SOURCES FXVTKCanvas.cpp vtkFXRenderWindowInteractor.cpp)

IF (WIN32)
	SET(VTKFOX_SOURCES ${VTKFOX_SOURCES} windows.cpp)
ELSE (WIN32)
	SET(VTKFOX_SOURCES ${VTKFOX_SOURCES} unix.cpp)
ENDIF (WIN32)

# suppress warnings
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -Wno-old-style-cast -Wno-deprecated -Wno-unused-parameter -Wno-long-long  -Wno-cast-qual -Wno-redundant-decls )
  REMOVE_DEFINITIONS ( -pedantic )
ENDIF ( )

ADD_LIBRARY(vtkfox ${VTKFOX_SOURCES})
