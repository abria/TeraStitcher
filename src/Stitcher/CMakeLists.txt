# CmakeLists.txt in Stitcher dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/IOManager)
include_directories (${TeraStitcher_SOURCE_DIR}/TinyXML)
include_directories (${TeraStitcher_SOURCE_DIR}/VolumeManager)
include_directories (${TeraStitcher_SOURCE_DIR}/CrossMIPs)


# Adds a library called stitcher (stitcher.a under Linux, stitcher.lib under Windows) from the all .h and .cpp files
file(GLOB stitcher_headers *.h)
file(GLOB stitcher_sources *.cpp)
add_library(stitcher STATIC ${stitcher_headers} ${stitcher_sources})