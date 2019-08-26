INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/src)

ADD_EXECUTABLE(cone cone.cpp)
TARGET_LINK_LIBRARIES(cone vtkfox ${FOX_LIBS} vtkRendering)

ADD_EXECUTABLE(xyplot xyplot.cpp)
TARGET_LINK_LIBRARIES(xyplot vtkfox ${FOX_LIBS} vtkHybrid)

ADD_EXECUTABLE(quad quad.cpp)
TARGET_LINK_LIBRARIES(quad vtkfox ${FOX_LIBS} vtkRendering)

ADD_EXECUTABLE(imgsurf imgsurf.cpp)
TARGET_LINK_LIBRARIES(imgsurf vtkfox ${FOX_LIBS} vtkRendering)

ADD_EXECUTABLE(simple simple.cpp)
TARGET_LINK_LIBRARIES(simple vtkfox ${FOX_LIBS} vtkRendering)

ADD_EXECUTABLE(double double.cpp)
TARGET_LINK_LIBRARIES(double vtkfox ${FOX_LIBS} vtkRendering)
