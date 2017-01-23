set(TeraStitcher_VERSION_MAJOR 1)
set(TeraStitcher_VERSION_MINOR 9)
set(TeraStitcher_VERSION_PATCH 66)
set(TeraStitcher_VERSION
  "${TeraStitcher_VERSION_MAJOR}.${TeraStitcher_VERSION_MINOR}.${TeraStitcher_VERSION_PATCH}")
add_definitions( -DTERASTITCHER_MAJOR=${TeraStitcher_VERSION_MAJOR} )
add_definitions( -DTERASTITCHER_MINOR=${TeraStitcher_VERSION_MINOR} )
add_definitions( -DTERASTITCHER_PATCH=${TeraStitcher_VERSION_PATCH} )