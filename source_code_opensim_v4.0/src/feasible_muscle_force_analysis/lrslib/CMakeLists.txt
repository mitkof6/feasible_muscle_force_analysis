# gmp
if (WIN32)
  include_directories(unistd_win/)
  include_directories(gmp_win64/)
  find_library(gmp
    NAMES libgmp-10
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/gmp_win64/
    PATH_SUFFIXES lib)
  set(DEPENDENCY_LIBRARIES ${gmp})
else ()
  set(DEPENDENCY_LIBRARIES gmp)
endif()

# lrslib
set(target lrs)
file(GLOB sources lrslib.c lrsgmp.c)
file(GLOB includes lrslib.h lrsgmp.h)
add_library(${target} STATIC ${sources} ${includes})
target_link_libraries (${target} ${DEPENDENCY_LIBRARIES})
target_compile_options(${target} PRIVATE -fPIC)
set_target_properties(
  ${target} PROPERTIES
  FOLDER "Libraries"
  PROJECT_LABEL "Libraries - ${target}")

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
