cmake_minimum_required(VERSION 3.1.0)

project(Spektrum)

set(CMAKE_AUTOMOC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(Spektrum STATIC
    include/exportlineprofiledialog.h
    include/feldialog.h
    include/linedialog.h
    include/lineprofilefitwindow.h
    include/linewindowbase.h
    include/selectlinetabledialog.h
    include/spectlist.h
    include/Spektrum.h
    exportlineprofiledialog.cpp
    feldialog.cpp
    linedialog.cpp
    lineprofilefitwindow.cpp
    linewindowbase.cpp
    selectlinetabledialog.cpp
    spectlist.cpp
    Spektrum.cpp
)

target_include_directories(Spektrum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Spektrum
    FCF
    FitObject
    Progressions
    SpectSimulation
    TermPlot
    Qt6::Widgets
)
