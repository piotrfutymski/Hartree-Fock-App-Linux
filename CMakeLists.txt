cmake_minimum_required(VERSION 3.10.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


SET ( CMAKE_CXX_FLAGS "-pthread" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_HOME_DIRECTORY}/build")
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_HOME_DIRECTORY}/build")
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_HOME_DIRECTORY}/build")


project(HartreeFockApp VERSION 0.1.0)

include(CTest)
enable_testing()

option(USE_CUDA "Use CUDA" ON)
include_directories("include" "CUDA")
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")
file(GLOB GPU_SOURCES "CUDA/*.cu")



if(USE_CUDA)
set(CMAKE_CUDA_COMPILER "/usr/local/cuda-11.1/bin/nvcc")
enable_language("CUDA")
    add_executable(HartreeFockApp ${SOURCES}  ${GPU_SOURCES} ${HEADERS})
else(USE_CUDA)
    add_executable(HartreeFockApp ${SOURCES} ${HEADERS})
endif()

message(STATUS "USE_CUDA: ${USE_CUDA}")
if(USE_CUDA)
    add_definitions(-DUSE_CUDA)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)



target_link_libraries(HartreeFockApp armadillo)
