cmake_minimum_required(VERSION 3.15)

project(mesh)

# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

#activation des options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -fPIC")
#compilation flags
set(debug 0)
if (debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g") # option debug
else(debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(debug)

#adding OpenMP library
find_package(OpenMP)
if(OPENMP_FOUND)
  #message("OpenMP found.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
  message("OpenMP not fbjhound.")
endif(OPENMP_FOUND)

### Recurse to the source code
add_subdirectory(src)

# install
install(
  TARGETS mesh
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib)
