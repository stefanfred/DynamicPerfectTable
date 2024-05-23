cmake_minimum_required(VERSION 3.28)
project(DynamicPerfectHashing)

set(CMAKE_CXX_STANDARD 17)

add_executable(DynamicPerfectHashing main.cpp
        DynPerfectHash.hpp)
