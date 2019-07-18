add_subdirectory(feasible_muscle_force_analysis)

file(GLOB applications *.cpp)

# dependencies
include_directories(include/)
set(DEPENDENCY_LIBRARIES ${OpenSim_LIBRARIES} FeasibleMuscleForceAnalysis)

# applications
foreach(application_program ${applications})
  # NAME_WE stands for "name without extension"
  get_filename_component(APPLICATION_NAME ${application_program} NAME_WE)
  add_executable(${APPLICATION_NAME} ${application_program})
  target_link_libraries(${APPLICATION_NAME} ${DEPENDENCY_LIBRARIES})
  set_target_properties(${APPLICATION_NAME}
    PROPERTIES
    PROJECT_LABEL "Application - ${APPLICATION_NAME}"
    FOLDER "Applications")
endforeach()
