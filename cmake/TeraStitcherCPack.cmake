# set CPack options
set (CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../DESCRIPTION.txt)
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt)
set (CPACK_PACKAGE_DESCRIPTION "A tool for fast automatic 3D-stitching of teravoxel-sized microscopy images")
set (CPACK_PACKAGE_VERSION_MAJOR ${TeraStitcher_VERSION_MAJOR} )
set (CPACK_PACKAGE_VERSION_MINOR ${TeraStitcher_VERSION_MINOR} )
set (CPACK_PACKAGE_VERSION_PATCH ${TeraStitcher_VERSION_PATCH} )
set (CPACK_PACKAGE_CONTACT "Alessandro Bria (a.bria@unicas.it)")
if(WITH_QT5)
	set(CPACK_PACKAGE_NAME "TeraStitcher-Qt5-standalone")
elseif(WITH_QT4)
	set(CPACK_PACKAGE_NAME "TeraStitcher-Qt4-standalone")
else()
	set(CPACK_PACKAGE_NAME "TeraStitcher-standalone")
endif()
set (CPACK_PACKAGING_INSTALL_PREFIX "/usr/local")
set (CPACK_BUNDLE_NAME TeraStitcher)
string (CONCAT DESKTOP_LINK_NAME 
	"TeraStitcher " 
	${TeraStitcher_VERSION_MAJOR} 
	"." 
	${TeraStitcher_VERSION_MINOR} 
	"." 
	${TeraStitcher_VERSION_PATCH})
set (CPACK_PACKAGE_EXECUTABLES "terastitcher-gui" "${DESKTOP_LINK_NAME}")
set (CPACK_CREATE_DESKTOP_LINKS "terastitcher-gui")

configure_file("${CMAKE_CURRENT_LIST_DIR}/TeraStitcherCPackOptions.cmake.in"
  "${TeraStitcher_BINARY_DIR}/TeraStitcherCPackOptions.cmake" @ONLY)
set(CPACK_PROJECT_CONFIG_FILE
  "${TeraStitcher_BINARY_DIR}/TeraStitcherCPackOptions.cmake")

include (CPack )