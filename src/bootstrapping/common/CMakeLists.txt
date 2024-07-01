message(STATUS "include common directory")

# Source files in this directory
set(COMMON_SOURCE_FILES ${COMMON_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/Choosemax.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MinicompFunc.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MinicompRemez.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Point.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PolyUpdate.cpp
	${CMAKE_CURRENT_LIST_DIR}/func.cpp
	${CMAKE_CURRENT_LIST_DIR}/Polynomial.cpp
	${CMAKE_CURRENT_LIST_DIR}/Remez.cpp
	${CMAKE_CURRENT_LIST_DIR}/RemezApp.cpp
)

set(COMMON_SOURCE_FILES ${COMMON_SOURCE_FILES} PARENT_SCOPE)
set(COMMON_HEADER_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)