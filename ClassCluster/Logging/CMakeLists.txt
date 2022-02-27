cmake_minimum_required(VERSION 3.1.0)

project(Logging)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Logging STATIC
    include/logfunction.h
    include/logger.h
    include/logmodel.h
    include/logwindow.h
    logfunction.cpp
    logger.cpp
    logmodel.cpp
    logwindow.cpp
)

target_include_directories(Logging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Logging PRIVATE
    MdiChild
    TableWindow
    Qt5::Widgets
)
