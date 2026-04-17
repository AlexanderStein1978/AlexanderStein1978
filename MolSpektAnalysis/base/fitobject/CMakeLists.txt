cmake_minimum_required(VERSION 3.5.0)

project(FitObject)

set(CMAKE_AUTOMOC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(FitObject STATIC
    include/fitobject.h
    include/gaussian.h
    include/gaussianwithsaturation.h
    include/lineprofile.h
    include/lorentzian.h
    include/lorentzianwithsaturation.h
    fitobject.cpp
    gaussian.cpp
    gaussianwithsaturation.cpp
	lineprofile.cpp
    lorentzian.cpp
    lorentzianwithsaturation.cpp
)

target_include_directories(FitObject PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FitObject
    Tools
<<<<<<< HEAD
    Qt6::Widgets
=======
    Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
