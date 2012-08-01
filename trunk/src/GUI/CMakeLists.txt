# CmakeLists.txt in GUI dir

# Make sure the compiler can find include files of other modules. 
include_directories (${TeraStitcher_SOURCE_DIR}/IOManager)
include_directories (${TeraStitcher_SOURCE_DIR}/TinyXML)
include_directories (${TeraStitcher_SOURCE_DIR}/VolumeManager)
include_directories (${TeraStitcher_SOURCE_DIR}/CrossMIPs)
include_directories (${TeraStitcher_SOURCE_DIR}/Stitcher)
include_directories (${TeraStitcher_SOURCE_DIR}/../include/tclap)

# Adds executable called terastitcher (terastitcher.a under Linux, terastitcher.exe under Windows) from the all .h and .cpp files
file(GLOB GUI_headers *.h)
file(GLOB GUI_sources *.cpp)
add_executable(terastitcher ${GUI_headers} ${GUI_sources})

# Link the executable to the other modules. 
target_link_libraries(terastitcher stitcher )
target_link_libraries(terastitcher volumemanager )
target_link_libraries(terastitcher crossmips )
target_link_libraries(terastitcher iomanager )
target_link_libraries(terastitcher tinyxml )
target_link_libraries(terastitcher ${OpenCV_LIBS})