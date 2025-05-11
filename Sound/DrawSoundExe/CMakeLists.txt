cmake_minimum_required(VERSION 3.1.0)

project(DrawSoundExe)

set(CMAKE_AUTOMOC ON)

add_executable(DrawSoundExe
    main.cpp
)

target_link_libraries(DrawSoundExe
    DrawSound
)
