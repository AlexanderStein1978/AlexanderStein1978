<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

cmake_policy(SET CMP0057 NEW)

project(ClassCluster_test)

add_definitions(-DDATA_DIRECTORY="${CMAKE_SOURCE_DIR}/ClassCluster/Data")

set(CMAKE_AUTOMOC ON)

include(GoogleTest)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_executable(ClassCluster_test
    test_Calculation.cpp
    main.cpp
)

target_link_libraries(ClassCluster_test
    Potential
    Simulation
    Templates
    TestHelper
    Tools
    gtest
<<<<<<< HEAD
    Qt6::Widgets
=======
    Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)

gtest_add_tests(TARGET ClassCluster_test)
