<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(Network)

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Network REQUIRED)
=======
find_package(Qt5 COMPONENTS Network REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

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
	Logging
	Simulation
<<<<<<< HEAD
	Qt6::Network
=======
	Qt5::Network
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
