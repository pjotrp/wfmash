cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(wflign VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE) # does not work

#set(CMAKE_BUILD_TYPE Release)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Release Debug Generic." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_BUILD_TYPE} MATCHES Release)
    #set(EXTRA_FLAGS "-Ofast -march=x86-64-v3 -flto -fno-fat-lto-objects")
    # set(EXTRA_FLAGS "-Ofast -march=x86-64-v3 -g")
    #set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    # Debug use the defaults from parent - so commenting out:
    # set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -g -fsanitize=address")
    # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -g -fsanitize=address")
else()
    # Use all standard-compliant optimizations - always add these:
    # set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
    # set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
endif ()

if(WFA_PNG_AND_TSV)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWFA_PNG_TSV_TIMING")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWFA_PNG_TSV_TIMING")
endif ()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#
# WFA2-lib Build
# use CMakelists to build rather than Makefile

#
# WFlign library
#
# target_link_libraries(wflign rt z)

set(wflign_SOURCE
  deps/atomic_image.cpp
  deps/lodepng/lodepng.cpp
  deps/atomic_image.cpp
  src/wflign.cpp
  src/wflign_alignment.cpp
  src/wflign_patch.cpp
  src/wflign_swizzle.cpp
  src/rkmh.cpp
  src/murmur3.cpp
)

set(wflign_INCLUDE
    deps/
    deps/lodepng
    deps/PicoSHA2
    deps/robin-hood-hashing
    ../../../deps/WFA2-lib
    ../../../deps/WFA2-lib/bindings/cpp
)

add_library(wflign_static STATIC ${wflign_SOURCE})
add_library(wflign SHARED ${wflign_SOURCE})
set_target_properties(wflign_static PROPERTIES OUTPUT_NAME wflign)
set_target_properties(wflign PROPERTIES OUTPUT_NAME wflign)

target_include_directories(wflign_static PRIVATE ${wflign_INCLUDE})
target_include_directories(wflign PRIVATE ${wflign_INCLUDE})
