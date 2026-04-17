<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(DrawSoundExe)

set(CMAKE_AUTOMOC ON)

add_executable(DrawSoundExe
    main.cpp
)

target_link_libraries(DrawSoundExe
    DrawSound
)
