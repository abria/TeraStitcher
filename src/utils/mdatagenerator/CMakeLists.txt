# CmakeLists.txt in mdatagenerator dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/imagemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/tclap)

link_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/libs)

# Adds executable called mdatagenerator (mdatagenerator under Linux, mdatagenerator.exe under Windows) from the all .h and .cpp files
file(GLOB mdatagenerator_headers *.h)
file(GLOB mdatagenerator_sources *.cpp)
add_executable(mdatagenerator ${mdatagenerator_headers} ${mdatagenerator_sources})

# Link the executable to the other modules. 
target_link_libraries(mdatagenerator imagemanager )
target_link_libraries(mdatagenerator iomanager )
target_link_libraries(mdatagenerator volumemanager )
target_link_libraries(mdatagenerator tinyxml )
target_link_libraries(mdatagenerator stitcher )
target_link_libraries(mdatagenerator crossmips )