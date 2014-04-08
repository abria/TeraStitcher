# CmakeLists.txt in VolumeManager dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/IOManager)
include_directories (${TeraStitcher_SOURCE_DIR}/Stitcher)
include_directories (${TeraStitcher_SOURCE_DIR}/TinyXML)
include_directories (${TeraStitcher_SOURCE_DIR}/../include)

# Adds a library called volumemanager (volumemanager.a under Linux, volumemanager.lib under Windows) from the all .h and .cpp files
file(GLOB volumemanager_headers *.h)
file(GLOB volumemanager_sources *.cpp)
add_library(volumemanager STATIC ${volumemanager_headers} ${volumemanager_sources})