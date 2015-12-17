# CMakeLists.txt in imagemanager dir

# add used modules to include path 
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)
include_directories (${TeraStitcher_SOURCE_DIR}/volumemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/tinyxml)

# add 3rd party modules to include path
include_directories (${CMAKE_SOURCE_DIR}/3rdparty/zlib)
include_directories (${CMAKE_BINARY_DIR}/3rdparty/zlib)
include_directories (${CMAKE_SOURCE_DIR}/3rdparty/libtiff)
include_directories (${CMAKE_BINARY_DIR}/3rdparty/libtiff)


# check HDF5 option
option (WITH_HDF5 "Support to HDF5" OFF) 
if(WITH_HDF5)
	message(STATUS "Checking HDF5")
	find_package (HDF5 NAMES hdf5 COMPONENTS C static)
	if(HDF5_FOUND)
		INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
		message(STATUS "Checking HDF5 ... OK!")
		add_definitions(-DENABLE_BDV_HDF5)
		# override hdf5 target properties to remove target zlib-target that causes conflicts with zlib already present in this project
		set_target_properties(${HDF5_C_STATIC_LIBRARY} PROPERTIES INTERFACE_LINK_LIBRARIES "szip-static")
		add_dependencies(${HDF5_C_STATIC_LIBRARY} zlib)
	else(HDF5_FOUND)
		message(FATAL_ERROR "Checking HDF5 ... error!!!")
		#message(FATAL_ERROR "Checking HDF5")
	endif(HDF5_FOUND)
endif(WITH_HDF5)

# add STATIC library from the all .h and .cpp files
file(GLOB imagemanager_headers *.h)
file(GLOB imagemanager_sources *.cpp)
add_library(imagemanager STATIC ${imagemanager_headers} ${imagemanager_sources})
target_link_libraries(imagemanager ${HDF5_C_STATIC_LIBRARY})

