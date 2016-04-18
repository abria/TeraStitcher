# set CPack options
set (CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../DESCRIPTION.txt)
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt)
set (CPACK_PACKAGE_DESCRIPTION "A tool for fast automatic 3D-stitching of teravoxel-sized microscopy images")
set (CPACK_PACKAGE_VERSION_MAJOR ${TeraStitcher_VERSION_MAJOR} )
set (CPACK_PACKAGE_VERSION_MINOR ${TeraStitcher_VERSION_MINOR} )
set (CPACK_PACKAGE_VERSION_PATCH ${TeraStitcher_VERSION_PATCH} )
set (CPACK_PACKAGE_CONTACT "Alessandro Bria (a.bria@unicas.it)")
set (CPACK_BUNDLE_NAME TeraStitcher)

configure_file("${CMAKE_CURRENT_LIST_DIR}/TeraStitcherCPackOptions.cmake.in"
  "${TeraStitcher_BINARY_DIR}/TeraStitcherCPackOptions.cmake" @ONLY)
set(CPACK_PROJECT_CONFIG_FILE
  "${TeraStitcher_BINARY_DIR}/TeraStitcherCPackOptions.cmake")

include (CPack )