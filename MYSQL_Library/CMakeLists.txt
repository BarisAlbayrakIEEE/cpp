# Set the minimum cmake version
cmake_minimum_required (VERSION 3.0)

# Set the project name
project (MYSQL_Library)

# Set the version number
set (MYSQL_Library_VERSION_MAJOR 1)
set (MYSQL_Library_VERSION_MINOR 0)

# Create a variable for the include directories
set (
	_INCLUDE_DIRECTORIES
	"${PROJECT_SOURCE_DIR}/source"
	"C:/Program Files/MySQL/MySQL Server 8.1/include"
	"C:/Program Files/MySQL/MySQL Connector C++ 8.0/include"
)

# Create a variable for the header files
set (
	_HEADER_FILES
	"${PROJECT_BINARY_DIR}/MYSQL_Library_Config.h"
	"${PROJECT_SOURCE_DIR}/source/headerMYSQL.hxx"
)	

# Create a variable for the source files
set (
	_SOURCE_FILES
	"${PROJECT_SOURCE_DIR}/source/main.cpp" 
	"${PROJECT_SOURCE_DIR}/source/sourceMYSQL.cxx" 
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
	"libmysql.lib"
	"mysqlcppconn8.lib"
)

# Create the config header file in the binary tree
configure_file (
	"${PROJECT_SOURCE_DIR}/MYSQL_Library_Config.h.in"
	"${PROJECT_BINARY_DIR}/MYSQL_Library_Config.h"
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
	MYSQL_Library 
	${_SOURCE_FILES}
)

# Create the executable
#add_executable(
#	MYSQL_Library 
#	${_SOURCE_FILES}
#)

# Link the libraries
target_link_libraries(
	MYSQL_Library
	${_STATIC_DEPENDENCIES}
)

# Set the output directories (DEBUG and RELEASE)
set_target_properties(
	MYSQL_Library
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
)
