# CmakeLists.txt in UI dir

# add used modules to include path
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)
include_directories (${TeraStitcher_SOURCE_DIR}/tinyxml)
include_directories (${TeraStitcher_SOURCE_DIR}/volumemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/crossmips)
include_directories (${TeraStitcher_SOURCE_DIR}/stitcher)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty)

# add executable from the all .h and .cpp files
file(GLOB example_headers *.h)
file(GLOB example_sources *.cpp)
add_executable(example_utility ${example_headers} ${example_sources})

# link the executable to the other modules 
target_link_libraries(example_utility volumemanager )
target_link_libraries(example_utility stitcher )
target_link_libraries(example_utility crossmips )
target_link_libraries(example_utility iomanager )
target_link_libraries(example_utility tinyxml )