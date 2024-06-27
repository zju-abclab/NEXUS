message(STATUS "include cnn directory")

# Source files in this directory
set(CNN_SOURCE_FILES ${CNN_SOURCE_FILES}
	${CMAKE_CURRENT_LIST_DIR}/cnn_seal.cpp
	${CMAKE_CURRENT_LIST_DIR}/infer_seal.cpp
)

set(CNN_SOURCE_FILES ${CNN_SOURCE_FILES} PARENT_SCOPE)
set(CNN_HEADER_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

