cmake_minimum_required(VERSION 3.1.0)

project(Network)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Network REQUIRED)

add_library(Network STATIC
	include/network.h
    include/networkclient.h
    include/networkserver.h
    network.cpp
    networkclient.cpp
	networkserver.cpp
)

target_include_directories(Network PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Network
	Simulation
	Qt5::Network
)
