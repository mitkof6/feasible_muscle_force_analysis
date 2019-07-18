# lrslib
if (WIN32)
  add_definitions(-DLRSLONG)
  file(GLOB lrs_sources
    lrslib/lrslib.c
    lrslib/lrslong.c)
  file(GLOB lrs_includes
    lrslib/lrslib.h
    lrslib/lrslong.h
    lrslib/unistd_win32/unistd.h)
else ()
  add_definitions(-DGMP)
  file(GLOB lrs_sources
    lrslib/lrslib.c
    lrslib/lrsgmp.c)
  file(GLOB lrs_includes
    lrslib/lrslib.h
    lrslib/lrsgmp.h)
  set(lrs_libraries gmp)
endif()
set(DEPENDENCY_LIBRARIES ${OpenSim_LIBRARIES} ${lrs_libraries})

# FeasibleMuscleForceAnalysis
file(GLOB sources *.cpp internal/*.cpp)
file(GLOB includes *.h internal/*.h)
set(target FeasibleMuscleForceAnalysis)
add_library(${target} SHARED
  ${sources} ${includes} ${lrs_sources} ${lrs_includes})
target_link_libraries (${target} ${DEPENDENCY_LIBRARIES})
set_target_properties(
  ${target} PROPERTIES
  FOLDER "Plugins"
  PROJECT_LABEL "Plugins - ${target}")

# install library
install(
  TARGETS ${target}
  EXPORT ${TARGET_EXPORT_NAME}
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

# install headers
set(INCLUDE_PREFIX "${CMAKE_INSTALL_INCLUDEDIR}")
install(FILES ${includes} DESTINATION ${INCLUDE_PREFIX}/${target})
