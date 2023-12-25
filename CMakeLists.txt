cmake_minimum_required(VERSION 3.10)
project(PrivateGPT VERSION 1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_BUILD_TYPE DEBUG)

find_package(SEAL 4.1)
add_executable(
    main 
    ${CMAKE_SOURCE_DIR}/src/main.cpp 
    ${CMAKE_SOURCE_DIR}/src/gelu.cpp 
    ${CMAKE_SOURCE_DIR}/src/layer_norm.cpp 
    ${CMAKE_SOURCE_DIR}/src/ckks_evaluator.cpp 
    ${CMAKE_SOURCE_DIR}/src/softmax.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_mul.cpp
)


target_link_libraries(main PRIVATE SEAL::seal)

# add_subdirectory(src)