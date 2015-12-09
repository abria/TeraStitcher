# CmakeLists.txt in UI dir

# add used modules to include path
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)
include_directories (${TeraStitcher_SOURCE_DIR}/tinyxml)
include_directories (${TeraStitcher_SOURCE_DIR}/volumemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/crossmips)
include_directories (${TeraStitcher_SOURCE_DIR}/stitcher)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty)
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/tclap)

link_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/libs)

# add executable from the all .h and .cpp files
file(GLOB UI_headers *.h)
file(GLOB UI_sources *.cpp)
add_executable(terastitcher ${UI_headers} ${UI_sources})

# link the executable to the other modules 
target_link_libraries(terastitcher stitcher )
target_link_libraries(terastitcher volumemanager )
target_link_libraries(terastitcher crossmips )
target_link_libraries(terastitcher iomanager )
target_link_libraries(terastitcher tinyxml )