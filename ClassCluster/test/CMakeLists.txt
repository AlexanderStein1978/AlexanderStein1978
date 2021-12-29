cmake_minimum_required(VERSION 3.1.0)

cmake_policy(SET CMP0057 NEW)

project(ClassCluster_test)

set(CMAKE_AUTOMOC ON)

include(GoogleTest)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(ClassCluster_test
    test_Calculation.cpp
    main.cpp
)

target_link_libraries(ClassCluster_test
    Potential
    Simulation
    Templates
    Tools
    gtest
    Qt5::Widgets
)

gtest_add_tests(TARGET ClassCluster_test)
