cmake_minimum_required(VERSION 3.11.0)
project(p2ptrans)

#add_library(lapjv source/lapjv.cpp)
#add_library(utils source/utils.f90)

find_package(F2PY REQUIRED)

set(name "fmodules")
set(f_sources "${CMAKE_CURRENT_SOURCE_DIR}/source/lap.f90 ${CMAKE_CURRENT_SOURCE_DIR}/source/tiling.f90 ${CMAKE_CURRENT_SOURCE_DIR}/source/potential.f90 ${CMAKE_CURRENT_SOURCE_DIR}/source/transform.f90 ${CMAKE_CURRENT_SOURCE_DIR}/source/lapjv.f90")
set(functions "munkres sphere circle distance derivative intoptimization fastoptimization lapjv")
set(flags --f90flags="-g -Og -fbacktrace -fopenmp" -lgomp)

set(generated_module_file "${CMAKE_CURRENT_BINARY_DIR}/${name}${PYTHON_EXTENSION_MODULE_SUFFIX}")

message(STATUS "foo include dir: ${generated_module_file}")

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
  )

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE}
    -m ${name}
    -c
 #   -l ${LIBRARY_OUTPUT_DIRECTORY}
    ${flags}
    ${f_sources}
    only: ${functions}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )



install(FILES bob DESTINATION p2ptrans)
