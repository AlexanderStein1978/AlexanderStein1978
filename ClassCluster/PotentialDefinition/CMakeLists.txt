<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(PotentialDefinition)

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets Network REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets Network REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_library(PotentialDefinition STATIC
    include/potentialdefinerinputdata.h
    include/potentialdefiner.h
    potentialdefinerinputdata.cpp
    potentialdefiner.cpp
)

target_include_directories(PotentialDefinition PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialDefinition
    Simulation
)
