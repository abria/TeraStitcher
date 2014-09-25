# Top-Level CmakeLists.txt to build TeraStitcher
project( TeraStitcher )
cmake_minimum_required(VERSION 2.7)

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


# disable annoying warnings on MSVC compilers
warnings_disable(CMAKE_CXX_FLAGS /wd4290) # vs2012: can't properly handle exceptions
warnings_disable(CMAKE_CXX_FLAGS /wd4996) # vs2012: complains about unsafe standard C++ functions

# build modules
add_subdirectory( 3rdparty/zlib )
add_subdirectory( 3rdparty/libtiff )
add_subdirectory( tinyxml )
add_subdirectory( imagemanager )
add_subdirectory( iomanager )
add_subdirectory( crossmips )
add_subdirectory( volumemanager )
add_subdirectory( stitcher )
add_subdirectory( ui )

# build user-selectable modules
# --- find modules ---
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/utils)
# --- ask for selection of the modules to be included in the project (default is OFF) ---
FOREACH(subdir ${SUBDIRS})
	option (WITH_UTILITY_MODULE_${subdir} "Include this plugin in the project" OFF) 
ENDFOREACH()
# --- include selected modules only ---
FOREACH(subdir ${SUBDIRS})
	if(WITH_UTILITY_MODULE_${subdir})
		message(STATUS "Include ./utils/" "${subdir} ...")
		add_subdirectory(utils/${subdir})
		message(STATUS "Include ./utils/" "${subdir} ... OK!")
	endif(WITH_UTILITY_MODULE_${subdir})
ENDFOREACH()
