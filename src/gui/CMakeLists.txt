# CMakeLists.txt in GUI dir

# add used modules to include path
include_directories (${TeraStitcher_SOURCE_DIR}/iomanager)
include_directories (${TeraStitcher_SOURCE_DIR}/imagemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/tinyxml)
include_directories (${TeraStitcher_SOURCE_DIR}/volumemanager)
include_directories (${TeraStitcher_SOURCE_DIR}/crossmips)
include_directories (${TeraStitcher_SOURCE_DIR}/stitcher)
include_directories (${TeraStitcher_SOURCE_DIR}/common)

# add 3rd party modules to include path
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty)
include_directories(${TeraStitcher_SOURCE_DIR}/3rdparty/tclap)


# add executable from the all .h and .cpp files
file(GLOB GUI_headers *.h)
file(GLOB GUI_sources *.cpp)

# set CMake 3.* automoc on, so that CMake automatically determines the source files that need moc pre-processing
set(CMAKE_AUTOMOC ON)

# set CMake 3.* autorcc on, so that CMake automatically determines the source files that need rcc pre-processing
set(CMAKE_AUTORCC ON)

# check for either Qt5 or Qt4
option (WITH_QT5 "With QT5" OFF) 
if(WITH_QT5)
	find_package(Qt5Widgets REQUIRED)
	add_executable(terastitcher-gui ${GUI_headers} ${GUI_sources} icons.qrc)
	target_link_libraries(terastitcher-gui Qt5::Widgets)
endif(WITH_QT5)
option (WITH_QT4 "With QT4" OFF) 
if(WITH_QT4)
	find_package(Qt4 REQUIRED)
	add_executable(terastitcher-gui ${GUI_headers} ${GUI_sources} icons.qrc)
	target_link_libraries(terastitcher-gui Qt4::QtGui)
endif(WITH_QT4)
if(NOT WITH_QT4 AND NOT WITH_QT5)
	message(FATAL_ERROR "Please select at least one Qt library (Qt4 or Qt5)")
endif(NOT WITH_QT4 AND NOT WITH_QT5)

# add executable to install target, including dependencies (Qt libraries)
SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/terastitcher-gui")
IF(WIN32)
  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/terastitcher-gui.exe")
ENDIF(WIN32)
INSTALL(TARGETS terastitcher-gui 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
    )
SET(DIRS ${QT_LIBRARY_DIRS})
INSTALL(CODE "
    include(BundleUtilities)
    fixup_bundle( \"${APPS}\" \"\" \"${DIRS}\" )
    " COMPONENT Runtime)

#target internal dependencies
target_link_libraries(terastitcher-gui imagemanager )
target_link_libraries(terastitcher-gui iomanager )
target_link_libraries(terastitcher-gui volumemanager )
target_link_libraries(terastitcher-gui tinyxml )
target_link_libraries(terastitcher-gui stitcher )
target_link_libraries(terastitcher-gui crossmips )
target_link_libraries(terastitcher-gui common )