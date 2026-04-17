<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(TestHelper)

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_library(TestHelper STATIC
    include/CalculationTestHelper.h
    CalculationTestHelper.cpp
)

target_include_directories(TestHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TestHelper
    Simulation
<<<<<<< HEAD
    Qt6::Widgets
=======
    Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
