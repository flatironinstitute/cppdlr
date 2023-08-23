# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# --------------------------------------------------------
# Build & Run the C++ doc examples and capture the output
# --------------------------------------------------------

add_custom_target(${PROJECT_NAME}_docs_example_output)
file(GLOB_RECURSE ExampleList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach(example ${ExampleList})
  get_filename_component(f ${example} NAME_WE)
  get_filename_component(d ${example} DIRECTORY)
  add_executable(${PROJECT_NAME}_doc_${f} EXCLUDE_FROM_ALL ${example})
  set_property(TARGET ${PROJECT_NAME}_doc_${f} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${d})
  add_custom_command(TARGET ${PROJECT_NAME}_doc_${f}
   COMMAND ${PROJECT_NAME}_doc_${f} > ${CMAKE_CURRENT_SOURCE_DIR}/${d}/${f}.output 2>/dev/null
   WORKING_DIRECTORY ${d}
  )
  add_dependencies(${PROJECT_NAME}_docs_example_output ${PROJECT_NAME}_doc_${f})
endforeach()

# ---------------------------------
# Top Sphinx target
# ---------------------------------
if(NOT DEFINED SPHINXBUILD_EXECUTABLE)
  find_package(Sphinx)
endif()

# Sphinx has internal caching, always run it
add_custom_target(${PROJECT_NAME}_docs_sphinx ALL)
add_custom_command(
  TARGET ${PROJECT_NAME}_docs_sphinx
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:${h5_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SPHINXBUILD_EXECUTABLE} -c . -j auto -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

option(Sphinx_Only "When building the documentation, skip the Python Modules and the generation of C++ Api and example outputs" OFF)
if(NOT Sphinx_Only)
  # Generation of C++ Api and Example Outputs
  add_dependencies(${PROJECT_NAME}_docs_sphinx ${PROJECT_NAME}_docs_example_output)
endif()


#----------------------------------
# Doxygen
#----------------------------------

find_package(Doxygen REQUIRED)

# Set input and output directories
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/c++)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)

# Configure Doxyfile
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

# Add custom target for Doxygen
add_custom_target(doxygen ALL
  #    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating Doxygen documentation"
  VERBATIM
)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/${PROJECT_NAME}
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|ico|js|xsl|css|py|txt|inv|bib|ttf|woff2|eot|sh)$"
  PATTERN "_*"
)
