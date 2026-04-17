<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(Logging)

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_library(Logging STATIC
    include/logfunction.h
    include/logger.h
    include/loglistitem.h
    include/loglist.h
    include/logmodel.h
    include/logwindow.h
    logfunction.cpp
    logger.cpp
    loglist.cpp
    logmodel.cpp
    logwindow.cpp
)

target_include_directories(Logging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Logging PRIVATE
    MdiChild
	Network
    TableWindow
<<<<<<< HEAD
    Qt6::Widgets
=======
    Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
