cmake_minimum_required(VERSION 3.0.0)
project(function-sampler VERSION 0.1.0)

add_library(function-sampler INTERFACE)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

cmake_minimum_required(VERSION 3.0)
project(function-sampler)

target_include_directories(function-sampler INTERFACE
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(function-sampler INTERFACE cxx_std_20)