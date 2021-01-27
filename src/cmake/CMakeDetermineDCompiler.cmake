# Find the D compiler
find_program(
    CMAKE_D_COMPILER
        NAMES "dmd"
        HINTS "${CMAKE_SOURCE_DIR}"
        DOC "DMD64 compiler"
)
mark_as_advanced(CMAKE_D_COMPILER)

set (CMAKE_D_SOURCE_FILE_EXTENSIONS d)
if (WIN32)
  set (CMAKE_D_OUTPUT_EXTENSION .obj)
else ()
  set (CMAKE_D_OUTPUT_EXTENSION .o)
endif ()
set (CMAKE_D_COMPILER_ENV_VAR "D")

# Configure variables set in this file for fast reload later on
configure_file(${CMAKE_CURRENT_LIST_DIR}/CMakeDCompiler.cmake.in
               ${CMAKE_PLATFORM_INFO_DIR}/CMakeDCompiler.cmake)
