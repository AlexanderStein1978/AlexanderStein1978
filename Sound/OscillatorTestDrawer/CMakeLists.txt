cmake_minimum_required(VERSION 3.1.0)

project(OscillatorTestDrawer)

set(CMAKE_AUTOMOC ON)

add_executable(OscillatorTestDrawer
    main.cpp
)

# target_include_directories(DrawSound PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(OscillatorTestDrawer
    DrawSound
)
