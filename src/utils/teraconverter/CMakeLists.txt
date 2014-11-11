# CmakeLists.txt in teraconverter dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/imagemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)
include_directories (${TeraStitcher_SOURCE_DIR}/utils/volumeconverter)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/tclap)

# Adds executable called teraconverter (teraconverter under Linux, teraconverter.exe under Windows) from the all .h and .cpp files
file(GLOB teraconverter_headers *.h)
file(GLOB teraconverter_sources *.cpp)
add_executable(teraconverter ${teraconverter_headers} ${teraconverter_sources})

# Link the executable to the other modules. 
target_link_libraries(teraconverter volumeconverter )
target_link_libraries(teraconverter imagemanager )
target_link_libraries(teraconverter iomanager )
