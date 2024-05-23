cmake_minimum_required(VERSION 3.28)
project(DynamicPerfectHashing)

set(CMAKE_CXX_STANDARD 20)

add_executable(DynamicPerfectHashing main.cpp
        DynPerfectHash.hpp)
