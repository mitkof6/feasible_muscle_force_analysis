CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

################################################################################

# project
project(FeasibleMuscleForceAnalysis VERSION 1.0.0)

# for rdm (emacs)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# for grouping into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# ctest
include(CTest)

if(UNIX AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the build type" FORCE)
endif()


# c++11, -g option is used to export debug symbols for gdb
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU OR
    ${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  # Using C++11 on OSX requires using libc++ instead of libstd++.
  # libc++ is an implementation of the C++ standard library for OSX.
  if(APPLE)
    if(XCODE)
      set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
      set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
    endif()
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
  endif()
endif()

# report platform
if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
  set(PLATFORM_ABI x64)
else()
  set(PLATFORM_ABI x86)
endif()
set(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" CACHE STRING
  "This is the platform and ABI we're building for." FORCE)

# install layout
# Layout. This works for all platforms:
#   * <prefix>/lib/cmake/<PROJECT-NAME>
#   * <prefix>/lib/
#   * <prefix>/include/
#   * <prefix>/bin/
#   * <prefix>/doc/
#   * <prefix>/scripts/
set(CMAKE_INSTALL_BINDIR bin)
set(CMAKE_INSTALL_INCLUDEDIR include)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_DOCDIR ${CMAKE_INSTALL_PREFIX}/docs)
set(CMAKE_INSTALL_CONFIG ${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME})

# define library postfix so that different builds will produce
# distinguished libraries
set(CMAKE_RELEASE_POSTFIX "_r" CACHE STRING "Release postfix")
set(CMAKE_DEBUG_POSTFIX "_d" CACHE STRING "Debug postfix")
set(CMAKE_RELWITHDEBINFO_POSTFIX "_rd" CACHE STRING
  "Release with debug info postfix")
set(CMAKE_MINSIZEREL_POSTFIX "_mr" CACHE STRING
  "Minimum size release postfix")
mark_as_advanced(CMAKE_RELEASE_POSTFIX)
mark_as_advanced(CMAKE_DEBUG_POSTFIX)
mark_as_advanced(CMAKE_RELWITHDEBINFO_POSTFIX)
mark_as_advanced(CMAKE_MINSIZEREL_POSTFIX)

# build executables and libraries in the same place
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR} CACHE INTERNAL "")
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR} CACHE INTERNAL "")

# configure project global settings that are include in the source
set(DATA_DIR "${PROJECT_SOURCE_DIR}/../data")
configure_file("${PROJECT_SOURCE_DIR}/src/Settings.h.in"
  "${PROJECT_SOURCE_DIR}/src/Settings.h")
set(CMAKE_EXTRA_INCLUDE_FILES)

################################################################################

# configure for new version of OpenSim (v4.0)
find_package(OpenSim REQUIRED)
if(OpenSim_FOUND)
    message(STATUS "OpenSim found at ${OpenSim_DIR}")
    include_directories(${OpenSim_INCLUDE_DIRS})
    link_directories(${OpenSim_LIB_DIR})
else()
    message(FATAL_ERROR "Can't find OpenSim at ${OpenSim_DIR}")
endif()

################################################################################

# submodules
add_subdirectory(src)
