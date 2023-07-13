# List of all examples
set(all_examples 
  #  myexample.cpp
  )

foreach(example ${all_examples})
  get_filename_component(example_name ${example} NAME_WE)
  get_filename_component(example_dir ${example} DIRECTORY)
  add_executable(${example_name} ${example})
  target_link_libraries(${example_name} ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings)
  set_property(TARGET ${example_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${example_dir})
  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${example_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()
  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${example_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++20
      --template=gcc
      --verbose
      --force
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${example}
    )
  endif()
endforeach()
