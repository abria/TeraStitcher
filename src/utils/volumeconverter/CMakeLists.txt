# CmakeLists.txt in volumeconverter dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/imagemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)

# add 3rd party modules to include path

# Adds a library called VolumeConverter (volumeconverter.a under Linux, volumeconverter.lib under Windows) from the all .h and .cpp files
file(GLOB volumeconverter_headers *.h)
file(GLOB volumeconverter_sources *.cpp)
add_library(volumeconverter STATIC ${volumeconverter_headers} ${volumeconverter_sources})
