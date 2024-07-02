cmake_minimum_required(VERSION 3.20)

project(NexusCUDA LANGUAGES CXX CUDA)

# CMake config
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

# CUDA config
if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
endif ()
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

# Set CUDA C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Print CUDA info
message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
message(STATUS "CMAKE_CUDA_COMPILER_ID = ${CMAKE_CUDA_COMPILER_ID}")
message(STATUS "CMAKE_CUDA_COMPILER_VERSION = ${CMAKE_CUDA_COMPILER_VERSION}")

# PhantomFHE config
option(PHANTOM_USE_CUDA_PTX "Use CUDA PTX Assembly" ON)
message(STATUS "Use CUDA PTX Assembly: ${PHANTOM_USE_CUDA_PTX}")
if (PHANTOM_USE_CUDA_PTX)
    add_compile_definitions(PHANTOM_USE_CUDA_PTX)
endif ()

# PhantomFHE
include_directories(thirdparty/phantom-fhe/include)

# NexusCUDA
file(GLOB_RECURSE NEXUS_CUDA_SOURCE_FILES CONFIGURE_DEPENDS
  "src/main.cu"
  "src/gelu.cu"
  "src/softmax.cu"
  "src/layer_norm.cu"
  # "src/matrix_mul.cu"
  "src/ckks_evaluator.cu"
  "thirdparty/phantom-fhe/src/*.cu"
)
add_executable(main ${NEXUS_CUDA_SOURCE_FILES})

# NexusCUDA Boostrapping
# include_directories(
#   /usr/local/include
#   /usr/local/include/NTL
#   src/bootstrapping/commmon
#   src/bootstrapping
# )

# file(GLOB_RECURSE NEXUS_CUDA_BOOTSTRAPING_SOURCE_FILES CONFIGURE_DEPENDS
#   "src/boostrap_test.cu"
#   "src/bootstrapping/*.cu"
#   "src/gelu.cu"
#   "src/softmax.cu"
#   "src/layer_norm.cu"
#   "src/matrix_mul.cu"
#   "src/ckks_evaluator.cu"
#   "thirdparty/phantom-fhe/src/*.cu"
# )

# add_executable(bootstrapping ${NEXUS_CUDA_BOOTSTRAPING_SOURCE_FILES})
