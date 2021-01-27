cmake_minimum_required (VERSION 3.0)
project (GDTK)
# CMakeLists.txt for the overall GasDynamic Toolkit
# We will do an out-of-source build by invoking from the build directory.
# cd dgd/build
# cmake ../src
# make
# make install

# message ("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
# message ("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
# find_path (MY_FOUND_PATH CMakeDetermineDCompiler.cmake ../src/cmake)
# message ("MY_FOUND_PATH: ${MY_FOUND_PATH}")
enable_language(D)

set (GAS_DIR ./gas)
include (${GAS_DIR}/gas_files.cmake)

set (ALL_FILES_D "")
foreach (item ${GAS_MODEL_FILES_D})
  set (ALL_FILES_D ${ALL_FILES_D} ${GAS_DIR}/${item})
endforeach ()

message ("Need to do some work...")
message ("ALL_FILES_D: ${ALL_FILES_D}")
