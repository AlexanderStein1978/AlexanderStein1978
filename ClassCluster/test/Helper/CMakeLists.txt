cmake_minimum_required(VERSION 3.1.0)

project(TestHelper)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TestHelper STATIC
    include/CalculationTestHelper.h
    CalculationTestHelper.cpp
)

target_include_directories(TestHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TestHelper
    Simulation
    Qt5::Widgets
)
