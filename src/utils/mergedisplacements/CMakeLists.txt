# CmakeLists.txt in mergedisplacements dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/imagemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)
include_directories (${TeraStitcher_SOURCE_DIR}/tinyxml)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/tclap)

# Adds executable called mergedisplacements (mergedisplacements under Linux, mergedisplacements.exe under Windows) from the all .h and .cpp files
file(GLOB mergedisplacements_headers *.h)
file(GLOB mergedisplacements_sources *.cpp)
add_executable(mergedisplacements ${mergedisplacements_headers} ${mergedisplacements_sources})

# Link the executable to the other modules. 
target_link_libraries(mergedisplacements iomanager )
target_link_libraries(mergedisplacements tinyxml )

