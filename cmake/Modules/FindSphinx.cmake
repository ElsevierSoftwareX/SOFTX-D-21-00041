find_program(SPHINX_BUILD_EXECUTABLE
  NAMES sphinx-build)

# does not seem to work
#if(SPHINX_BUILD_EXECUTABLE)
#  execute_process(
#    COMMAND ${SPHINX_BUILD_EXECUTABLE} --version
#    OUTPUT_VARIABLE SPHINX_VERSION
#    )
#  string(REPLACE "sphinx-build" "" SPHINX_VERSION "${SPHINX_VERSION}")
#  string(STRIP "${SPHINX_VERSION}" SPHINX_VERSION)
#endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sphinx
  REQUIRED_VARS SPHINX_BUILD_EXECUTABLE
#  VERSION_VAR SPHINX_VERSION
)

mark_as_advanced(SPHINX_BUILD_EXECUTABLE)
