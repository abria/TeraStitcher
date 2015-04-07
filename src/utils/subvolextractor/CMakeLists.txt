# CmakeLists.txt in subvolextractor dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/imagemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/tclap)

# Adds executable called subvolextractor (subvolextractor under Linux, subvolextractor.exe under Windows) from the all .h and .cpp files
file(GLOB subvolextractor_headers *.h)
file(GLOB subvolextractor_sources *.cpp)
add_executable(subvolextractor ${subvolextractor_headers} ${subvolextractor_sources})

# Link the executable to the other modules. 
target_link_libraries(subvolextractor imagemanager )
target_link_libraries(subvolextractor iomanager )
