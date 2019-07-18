# dependencies
if (WIN32)
    include_directories(unistd_win32/)
    # include_directories(gmp_win32/)
    #link_directories(${CMAKE_CURRENT_SOURCE_DIR}/gmp_win32/)
    #find_library(gmp
    #             NAMES gmp
    #             PATHS ${CMAKE_CURRENT_SOURCE_DIR}/gmp_win32/
    #             PATH_SUFFIXES lib)
endif()
# add_compile_definitions(LRSLONG)
# add_compile_definitions(GMP) # may not work with old versions
# add_definitions(-DGMP)
add_definitions(-DLRSLONG)
include_directories(lrslib/)
set(DEPENDENCY_LIBRARIES ${OpenSim_LIBRARIES} ${gmp})

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

# install headers
set(INCLUDE_PREFIX "${CMAKE_INSTALL_INCLUDEDIR}")
install(FILES ${includes} DESTINATION ${INCLUDE_PREFIX}/${target})
