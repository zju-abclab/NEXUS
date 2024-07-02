cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
project(PrivateGPT VERSION 1.0)

find_package(SEAL 4.1 REQUIRED)
add_executable(
    main
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/gelu.cpp
    ${CMAKE_SOURCE_DIR}/src/layer_norm.cpp
    ${CMAKE_SOURCE_DIR}/src/ckks_evaluator.cpp
    ${CMAKE_SOURCE_DIR}/src/softmax.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_mul.cpp
)

add_subdirectory(src/bootstrapping/common)
add_subdirectory(src/bootstrapping)

add_executable(
    bootstrapping
    ${CMAKE_SOURCE_DIR}/src/bootstrap_test.cpp
    ${CMAKE_SOURCE_DIR}/src/gelu.cpp
    ${CMAKE_SOURCE_DIR}/src/layer_norm.cpp
    ${CMAKE_SOURCE_DIR}/src/ckks_evaluator.cpp
    ${CMAKE_SOURCE_DIR}/src/softmax.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_mul.cpp
    ${COMMON_SOURCE_FILES}
    ${BOOTSTRAPPING_SOURCE_FILES}
)

target_include_directories(bootstrapping PUBLIC /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${BOOTSTRAPPING_HEADER_DIR})

target_link_libraries(main PRIVATE SEAL::seal)
target_link_libraries(bootstrapping PRIVATE ntl gmp m pthread SEAL::seal)

# add_subdirectory(src)
