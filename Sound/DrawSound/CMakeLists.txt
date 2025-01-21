cmake_minimum_required(VERSION 3.1.0)

project(DrawSound)

add_definitions(-DDATA_DIRECTORY="${CMAKE_SOURCE_DIR}/Sound/Data")

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets Multimedia REQUIRED)

add_executable(DrawSound
    main.cpp
    boxfilterdialog.cpp
    boxfilterdialog.h
    frequencywindow.cpp
    frequencywindow.h
    maxbuffer.cpp
    maxbuffer.h
    recordanddrawControl.cpp
    recordanddrawControl.h
    roundbuffer.h
    soundmainwindow.cpp
    soundmainwindow.h
    windowselectdialog.cpp
    windowselectdialog.h
    sounddrawwindow.cpp
    sounddrawwindow.h
    soundwindow.cpp
    soundwindow.h
    soundvector.h
    soundvector.cpp
    soundmatrix.h
    soundmatrix.cpp
)

# target_include_directories(DrawSound PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DrawSound
    DiagWindow
    Qt5::Widgets
    Qt5::Multimedia
)
