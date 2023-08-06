# Set the minimum cmake version
cmake_minimum_required (VERSION 3.0)

# Set the project name
project (GeometryLibrary_OCCT)

# Set the version number
set (GeometryLibrary_OCCT_VERSION_MAJOR 1)
set (GeometryLibrary_OCCT_VERSION_MINOR 0)

# Create a variable for the include directories
set (
	_INCLUDE_DIRECTORIES
	"${PROJECT_SOURCE_DIR}/source"
	"C:\\OpenCASCADE-7.6.0-vc14-64\\opencascade-7.6.0\\inc"
)

# Create a variable for the header files
set (
	_HEADER_FILES
	"${PROJECT_BINARY_DIR}/GeometryLibrary_OCCT_Config.h"
	"${PROJECT_SOURCE_DIR}/source/Macros.h"
	"${PROJECT_SOURCE_DIR}/source/GeometryAbstractObject.hxx"
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

# Create a variable for the source files
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

# Create a variable for the library directories
set (
	_LIBRARY_DIRECTORIES 
	"C:/OpenCASCADE-7.6.0-vc14-64/opencascade-7.6.0/win64/vc14/lib"
	"${PROJECT_SOURCE_DIR}/SharedLibraries"
)

# Create a variable for the static dependencies
set (
	_STATIC_DEPENDENCIES 
	"TKTopAlgo.lib"
	"TKernel.lib"
	"TKPrim.lib"
)

# Create the config header file in the binary tree
configure_file (
	"${PROJECT_SOURCE_DIR}/GeometryLibrary_OCCT_Config.h.in"
	"${PROJECT_BINARY_DIR}/GeometryLibrary_OCCT_Config.h"
)

# Add the include directories
# Add also PROJECT_BINARY_DIR to the include directories for Config.h
include_directories (
	${_INCLUDE_DIRECTORIES}
	"${PROJECT_BINARY_DIR}"
)

# Add the library directories
link_directories (
	${_LIBRARY_DIRECTORIES}
)

# Create the static library
add_library (
	GeometryLibrary_OCCT 
	${_SOURCE_FILES}
)

# Create the executable
#add_executable(
#	GeometryLibrary_OCCT 
#	${_SOURCE_FILES}
#)

# Link the libraries
target_link_libraries(
	GeometryLibrary_OCCT
	${_STATIC_DEPENDENCIES}
)

# Set the output directories (DEBUG and RELEASE)
set_target_properties(
	GeometryLibrary_OCCT
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)
