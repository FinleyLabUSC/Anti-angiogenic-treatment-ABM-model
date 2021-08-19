cmake_minimum_required(VERSION 3.0)
project(Angio_Treatment)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-fopenmp -std=c++14 -O3")
include_directories(inc)


file(GLOB SOURCES "src/*.cpp")

add_executable(runModel ${SOURCES})
