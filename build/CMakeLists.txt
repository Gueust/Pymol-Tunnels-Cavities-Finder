# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( TunnelsCavitiesFinder )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" 
  VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include_directories (BEFORE "../CGAL-Fix")
  include( CGAL_CreateSingleSourceCGALProgram )

   
  #include_directories (BEFORE	"/usr/local/include" 
	#			"/opt/ESBTL-1.0-beta01/include"
	#			"/usr/local/lib")
  #file(GLOB sources "../src/*.cpp" "../src/*.ui")
  #file(GLOB headers "../src/*.h")
  #add_executable(TunnelsCavitiesFinder ${sources} ${headers})

  create_single_source_cgal_program( "../src/TunnelsCavitiesFinder.cpp")

else()
  
    message(STATUS "This program requires the CGAL library, 
    and will not be compiled.")
  
endif()
