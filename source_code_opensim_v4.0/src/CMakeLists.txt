add_subdirectory(feasible_muscle_force_analysis)

file(GLOB tests test*.cpp)

# dependencies
include_directories(include/)
set(DEPENDENCY_LIBRARIES ${OpenSim_LIBRARIES} FeasibleMuscleForceAnalysis)

# tests
if(BUILD_TESTING)
  enable_testing()
  foreach(test_program ${tests})
    # NAME_WE stands for "name without extension"
    get_filename_component(TEST_NAME ${test_program} NAME_WE)
    add_executable(${TEST_NAME} ${test_program})
    target_link_libraries(${TEST_NAME}
      ${target} ${DEPENDENCY_LIBRARIES})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    set_target_properties(${TEST_NAME}
      PROPERTIES
      PROJECT_LABEL "Test - ${TEST_NAME}"
      FOLDER "Tests")
  endforeach()
endif()
