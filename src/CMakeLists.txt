cmake_minimum_required(VERSION 3.27)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 20)
set(COMPILER_FLAGS "-Wall --std=c++20 -O3 -march=native -funroll-loops")
set(LINKER_FLAGS "-static-libstdc++ -static-libgcc -static")
set(CMAKE_CXX_FLAGS "${COMPILER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LINKER_FLAGS}")
#set(CMAKE_CXX_CLANG_TIDY clang-tidy-12 -checks='*')

add_subdirectory(3rd)
include_directories(${CMAKE_CURRENT_LIST_DIR}/3rd/or-tools)
include_directories(${CMAKE_CURRENT_LIST_DIR}/3rd/cxxopts/include)

add_executable(mac_universal mac_universal.cpp ortools_kit.hpp aes.hpp sbox_espresso.hpp)
target_link_libraries(mac_universal ortools::ortools)

add_executable(smac_demo smac_demo.cpp smac_refcode.h smac_testvectors.h)
target_link_libraries(smac_demo)
