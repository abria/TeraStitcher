# find all headers and sources
file(GLOB headers *.h)
file(GLOB sources *.cpp)

# add parent dir to include path
include_directories (../../)

# get current directory name
get_filename_component(dirname ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# add static library
add_library(ioplugin_${dirname} STATIC ${headers} ${sources})
set_target_properties(ioplugin_${dirname} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# set up JNI library
find_package(JNI)
if (JNI_FOUND)
    message (STATUS "JNI_INCLUDE_DIRS=${JNI_INCLUDE_DIRS}")
    message (STATUS "JNI_LIBRARIES=${JNI_LIBRARIES}")
    include_directories (${JNI_INCLUDE_DIRS})
endif()

# link 3rd party libraries
target_link_libraries(ioplugin_${dirname} ${JNI_LIBRARIES})
