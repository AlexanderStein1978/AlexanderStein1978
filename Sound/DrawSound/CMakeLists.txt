cmake_minimum_required(VERSION 3.1.0)

project(DrawSound)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets Multimedia REQUIRED)

add_executable(DrawSound
    main.cpp
    recordanddrawControl.cpp
    recordanddrawControl.h
    sounddrawwindow.cpp
    sounddrawwindow.h
)

# target_include_directories(DrawSound PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DrawSound
    DiagWindow
    Qt5::Widgets
    Qt5::Multimedia
)
