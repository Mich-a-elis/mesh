cmake_minimum_required(VERSION 3.15)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

SET(SRCS
    mesh.cpp
    meshio.cpp
    vect.cpp
)

SET(INCLUDE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../include/")
SET(HEADERS
  ${INCLUDE_ROOT}/mesh.h
  ${INCLUDE_ROOT}/vect.h
)

add_library(mesh ${SRCS} ${HEADERS})

target_include_directories(mesh PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
