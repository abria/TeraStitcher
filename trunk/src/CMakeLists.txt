# Top-Level CmakeLists.txt to build TeraStitcher
project( TeraStitcher )
cmake_minimum_required(VERSION 2.7)

# hiding unused predefined variables
set (CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE INTERNAL "")

# setting up library and binary directories
if(NOT LIBRARY_OUTPUT_PATH)
  set (LIBRARY_OUTPUT_PATH ${TeraStitcher_BINARY_DIR}/libs CACHE INTERNAL "Single output directory for building all libraries.")
endif(NOT LIBRARY_OUTPUT_PATH)
if(NOT EXECUTABLE_OUTPUT_PATH)
  set (EXECUTABLE_OUTPUT_PATH ${TeraStitcher_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
endif(NOT EXECUTABLE_OUTPUT_PATH)

# setting up OpenCV library
find_package(OpenCV REQUIRED)
if(NOT OpenCV_FOUND)
  message(FATAL_ERROR "OpenCV library not found or not properly installed")
endif()
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIB_DIR})

# setting up "Release" as default build mode
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING 
        "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." 
        FORCE) # Default to Release
endif()

# build modules
add_subdirectory( TinyXML )
add_subdirectory( IOManager )
add_subdirectory( CrossMIPs )
add_subdirectory( VolumeManager )
add_subdirectory( Stitcher )
add_subdirectory( GUI )
