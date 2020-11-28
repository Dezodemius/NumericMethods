cmake_minimum_required(VERSION 3.17)
project(NumericMethods)

set(CMAKE_CXX_STANDARD 14)

add_executable(NumericMethods src/main.cpp)

file(MAKE_DIRECTORY pics)
file(MAKE_DIRECTORY gifs)
file(MAKE_DIRECTORY data)
file(MAKE_DIRECTORY bin)