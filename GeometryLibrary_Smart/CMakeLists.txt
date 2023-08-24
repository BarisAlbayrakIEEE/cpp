cmake_minimum_required (VERSION 3.0)
project (GeometryLibrary_Smart)
set (GeometryLibrary_Smart_VERSION_MAJOR 1)
set (GeometryLibrary_Smart_VERSION_MINOR 0)

# Include directories
set (
	_INCLUDE_DIRECTORIES
	"${PROJECT_SOURCE_DIR}/source"
)

# Library directories
#set (
#	_LIBRARY_DIRECTORIES 
#	"${PROJECT_SOURCE_DIR}/SharedLibraries"
#)

# Header files
set (
	_HEADER_FILES
	"${PROJECT_BINARY_DIR}/GeometryLibrary_Smart_Config.h"
	"${PROJECT_SOURCE_DIR}/source/Macros.h"
	"${PROJECT_SOURCE_DIR}/source/GeometryObject.hxx"
	"${PROJECT_SOURCE_DIR}/source/GeometryParameters.hxx"
	"${PROJECT_SOURCE_DIR}/source/GeometryException.hxx"
	"${PROJECT_SOURCE_DIR}/source/GeometryMath.hxx"
	"${PROJECT_SOURCE_DIR}/source/ReferenceAbstractObject.hxx"
	"${PROJECT_SOURCE_DIR}/source/ReferenceObject.hxx"
	"${PROJECT_SOURCE_DIR}/source/Axis.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Circle.hxx" 
	"${PROJECT_SOURCE_DIR}/source/CoordSystem.hxx" 
	"${PROJECT_SOURCE_DIR}/source/GlobalCoordSystem.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Line.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Plane.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Point2D.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Point3D.hxx" 
	"${PROJECT_SOURCE_DIR}/source/PointBase.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Vector2D.hxx" 
	"${PROJECT_SOURCE_DIR}/source/Vector3D.hxx" 
	"${PROJECT_SOURCE_DIR}/source/VectorBase.hxx"
)	

# Source files
set (
	_SOURCE_FILES
	"${PROJECT_SOURCE_DIR}/source/Axis.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Circle.cxx" 
	"${PROJECT_SOURCE_DIR}/source/CoordSystem.cxx" 
	"${PROJECT_SOURCE_DIR}/source/GeometryMath.cxx" 
	"${PROJECT_SOURCE_DIR}/source/GeometryObject.cxx" 
	"${PROJECT_SOURCE_DIR}/source/GlobalCoordSystem.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Line.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Plane.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Point2D.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Point3D.cxx" 
	"${PROJECT_SOURCE_DIR}/source/PointBase.cxx" 
	"${PROJECT_SOURCE_DIR}/source/ReferenceObject.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Vector2D.cxx" 
	"${PROJECT_SOURCE_DIR}/source/Vector3D.cxx" 
	"${PROJECT_SOURCE_DIR}/source/VectorBase.cxx"
)

# Static dependencies
#set (
#	_STATIC_DEPENDENCIES 
#	"TKTopAlgo.lib"
#	"TKernel.lib"
#	"TKPrim.lib"
#)

# Create the config header file in the binary tree
configure_file (
	"${PROJECT_SOURCE_DIR}/GeometryLibrary_Smart_Config.h.in"
	"${PROJECT_BINARY_DIR}/GeometryLibrary_Smart_Config.h"
)

# Add the include directories
# Add also PROJECT_BINARY_DIR to the include directories for Config.h
include_directories (
	${_INCLUDE_DIRECTORIES}
	"${PROJECT_BINARY_DIR}"
)

# Add the library directories
#link_directories (
#	${_LIBRARY_DIRECTORIES}
#)

# Create the static library
add_library (
	GeometryLibrary_Smart 
	${_SOURCE_FILES}
)

# Create the executable
#add_executable(
#	GeometryLibrary_Smart 
#	${_SOURCE_FILES}
#)

# Link the libraries
#target_link_libraries(
#	GeometryLibrary_Smart
#	${_STATIC_DEPENDENCIES}
#)

# Set the output directories (DEBUG and RELEASE)
set_target_properties(
	GeometryLibrary_Smart
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)
