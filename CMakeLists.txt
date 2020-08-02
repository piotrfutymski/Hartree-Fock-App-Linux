cmake_minimum_required(VERSION 3.10.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

SET ( CMAKE_CXX_FLAGS "-pthread" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_HOME_DIRECTORY}/build")
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_HOME_DIRECTORY}/build")
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_HOME_DIRECTORY}/build")


project(HartreeFockApp VERSION 0.1.0)

include_directories("include")
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")

include(CTest)
enable_testing()

add_executable(HartreeFockApp ${SOURCES} ${HEADERS})

target_link_libraries(HartreeFockApp armadillo)
