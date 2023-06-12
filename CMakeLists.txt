cmake_minimum_required(VERSION 3.10.2)

project(AirBreathing VERSION 1.1)

set(CMAKE_CXX_STANDARD 14)

file(GLOB_RECURSE SRC_LIB private/*.cpp)

include_directories(
	${CMAKE_SOURCE_DIR}/public
)

add_subdirectory(extern/pybind11)
pybind11_add_module(komurasakiairbreathing pythonmodule.cpp
	${SRC_LIB}
)
