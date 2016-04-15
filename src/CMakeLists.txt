# Top-Level CmakeLists.txt to build TeraStitcher
project( TeraStitcher )
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

# include TeraStitcher software version
include(version.cmake)

# include cmake scripts
include(cmake/utils.cmake)

# hide unused predefined variables
set (CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE INTERNAL "")

# set up library and binary directories
if(NOT LIBRARY_OUTPUT_PATH)
  set (LIBRARY_OUTPUT_PATH ${TeraStitcher_BINARY_DIR}/libs CACHE INTERNAL "Single output directory for building all libraries.")
endif(NOT LIBRARY_OUTPUT_PATH)
if(NOT EXECUTABLE_OUTPUT_PATH)
  set (EXECUTABLE_OUTPUT_PATH ${TeraStitcher_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
endif(NOT EXECUTABLE_OUTPUT_PATH)


# backward compatibility for MacOS X 10.5+
#SET(MAC_COVERAGE_COMPILE_FLAGS "-mmacosx-version-min=10.4")
#SET(MAC_COVERAGE_LINK_FLAGS    "-mmacosx-version-min=10.4")
#if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MAC_COVERAGE_COMPILE_FLAGS}" )
#	SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${MAC_COVERAGE_LINK_FLAGS}" )
#endif()

# disable annoying warnings on MSVC compilers
if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
	warnings_disable(CMAKE_CXX_FLAGS /wd4290) # vs2012: can't properly handle exceptions
	warnings_disable(CMAKE_CXX_FLAGS /wd4996) # vs2012: complains about unsafe standard C++ functions
endif()


# check GUI option
option (WITH_GUI "Build Qt-based GUI" OFF) 
if(WITH_GUI)
	message(STATUS "Include ./gui ...")
	add_subdirectory( gui )	
	message(STATUS "Include ./gui ... OK!")
endif(WITH_GUI)

# build modules
add_subdirectory( 3rdparty/zlib )
add_subdirectory( 3rdparty/libtiff )
add_subdirectory( tinyxml )
add_subdirectory( imagemanager )
add_subdirectory( iomanager )
add_subdirectory( crossmips )
add_subdirectory( volumemanager )
add_subdirectory( stitcher )
add_subdirectory( common )
add_subdirectory( ui )

# build user-selectable modules
# --- find modules ---
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/utils)
# --- ask for selection of the modules to be included in the project (default is OFF) ---
FOREACH(subdir ${SUBDIRS})
	option (WITH_UTILITY_MODULE_${subdir} "Include this utility module (optional) in the project" OFF) 
ENDFOREACH()
# --- include selected modules only ---
FOREACH(subdir ${SUBDIRS})
	if(WITH_UTILITY_MODULE_${subdir})
		message(STATUS "Include ./utils/" "${subdir} ...")
		add_subdirectory(utils/${subdir})
		message(STATUS "Include ./utils/" "${subdir} ... OK!")
	endif(WITH_UTILITY_MODULE_${subdir})
ENDFOREACH()

# set CPack
include (${CMAKE_CURRENT_SOURCE_DIR}/cpack/options.cmake)
include (CPack )
