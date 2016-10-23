cmake_minimum_required(VERSION 2.8)
project(posetofsets)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories("./"  ${Boost_INCLUDE_DIRS} )

add_executable(test test.cpp )


set(CMAKE_CXX_FLAGS "-O2   -fdiagnostics-color=auto -std=c++11 -m64  -Wall -Wextra -Wvla  " )

