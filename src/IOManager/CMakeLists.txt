# CmakeLists.txt in IOManager dir

# Adds a library called iomanager (iomanager.a under Linux, iomanager.lib under Windows) from the all .h and .cpp files
file(GLOB iomanager_headers *.h)
file(GLOB iomanager_sources *.cpp)
add_library(iomanager STATIC ${iomanager_headers} ${iomanager_sources})