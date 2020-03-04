CMAKE_MINIMUM_REQUIRED(VERSION 3.7)

PROJECT(BilateralFilter CXX)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)


############################## 3RD PARTIES LIBS ################################

FIND_PACKAGE(cpuid REQUIRED)
find_package( OpenCV HINTS "~/Software/OpenCVLibs/opencv-3.1.0/release/share/OpenCV/" )

INCLUDE_DIRECTORIES(
    SYSTEM
    ${CPUID_INCLUDE_DIRS}
    ${LIBS_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

LINK_LIBRARIES(
    ${CPUID_LIBRARIES}
    ${OpenCV_LIBS}
)


################################# PROJECT ##################################

SET(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/)

FILE(GLOB_RECURSE SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)
FILE(GLOB_RECURSE HEADER_FILES ${CMAKE_SOURCE_DIR}/include/*.hpp)

INCLUDE_DIRECTORIES(
    ${INCLUDE_DIR}
)

SET(TARGET_NAME bilteralFilter)
ADD_EXECUTABLE(${TARGET_NAME} ${SOURCE_FILES} ${HEADER_FILES})


TARGET_COMPILE_OPTIONS(
    ${TARGET_NAME} PRIVATE
    -g
    -O3
    -std=c++17
)

# Adding as much warning as possible on GNU gcc/g++ and Clang
IF (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    TARGET_COMPILE_OPTIONS(
        ${TARGET_NAME} PRIVATE
        -Wall
        -Wextra
    )
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "[Cc]lang")
    TARGET_COMPILE_OPTIONS(
        ${TARGET_NAME} PRIVATE
        -Weverything
        -Wfatal-errors
    )
    IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0.0")
    ELSE ()
    ENDIF ()
ENDIF ()


# Disable some warnings
TARGET_COMPILE_OPTIONS(
    ${TARGET_NAME} PRIVATE
    -Wno-unused-parameter
    -Wno-shadow
)
IF (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")

ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "[Cc]lang")
    TARGET_COMPILE_OPTIONS(
        ${TARGET_NAME} PRIVATE
        -Wno-padded
        -Wno-c++98-compat
        -Wno-weak-vtables
        -Wno-exit-time-destructors
    )
ENDIF ()
