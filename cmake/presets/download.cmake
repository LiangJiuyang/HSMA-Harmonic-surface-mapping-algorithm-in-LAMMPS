# preset that turns on packages with automatic downloads of sources of potentials
# compilation of libraries like Plumed or ScaFaCoS can take a considerable amount of time.

set(ALL_PACKAGES KIM LATTE MSCG VORONOI USER-PLUMED USER-SCAFACOS USER-SMD USER-MESONT)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(DOWNLOAD_KIM ON CACHE BOOL "" FORCE)
set(DOWNLOAD_LATTE ON CACHE BOOL "" FORCE)
set(DOWNLOAD_MSCG ON CACHE BOOL "" FORCE)
set(DOWNLOAD_VORO ON CACHE BOOL "" FORCE)
set(DOWNLOAD_EIGEN3 ON CACHE BOOL "" FORCE)
set(DOWNLOAD_PLUMED ON CACHE BOOL "" FORCE)
set(DOWNLOAD_SCAFACOS ON CACHE BOOL "" FORCE)

