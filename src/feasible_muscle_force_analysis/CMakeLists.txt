# find gmp
add_definitions(-DGMP)
if (WIN32)
  error("GMP library does not working properly on Windows")
else()
  find_path(GMP_INCLUDE_DIRS NAMES gmp.h)
  find_library(GMP_LIBRARY NAMES gmp libgmp)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GMP DEFAULT_MSG
    GMP_INCLUDE_DIRS
    GMP_LIBRARY)
  mark_as_advanced(GMP_INCLUDE_DIRS GMP_LIBRARY)
  include_directories(${GMP_INCLUDE_DIRS})
endif()
set(DEPENDENCY_LIBRARIES ${OpenSim_LIBRARIES} ${GMP_LIBRARY})

# FeasibleMuscleForceAnalysis
file(GLOB sources *.cpp internal/*.cpp lrslib/*.c)
file(GLOB includes *.h internal/*.h lrslib/*.h)
set(target FeasibleMuscleForceAnalysis)
add_library(${target} SHARED ${sources} ${includes})
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

# # install headers
# set(INCLUDE_PREFIX "${CMAKE_INSTALL_INCLUDEDIR}")
# install(FILES ${includes} DESTINATION ${INCLUDE_PREFIX}/${target})
