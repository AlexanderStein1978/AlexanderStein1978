<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(DrawSound)

add_definitions(-DDATA_DIRECTORY="${CMAKE_SOURCE_DIR}/Sound/Data")

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets Multimedia REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets Multimedia REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_library(DrawSound STATIC
    include/oscillator.h
    include/recordanddrawControl.h
    include/soundmainwindow.h
    boxfilterdialog.cpp
    boxfilterdialog.h
    definesamplesizedialog.cpp
    definesamplesizedialog.h
    frequencywindow.cpp
    frequencywindow.h
    maxbuffer.cpp
    maxbuffer.h
    patternbuffer.cpp
    patternbuffer.h
    recordanddrawControl.cpp
    roundbuffer.cpp
    roundbuffer.h
    soundmainwindow.cpp
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
    oscillator.cpp
    oscillatorarray.h
    oscillatorarray.cpp
    oscillatordiagram.h
    oscillatordiagram.cpp
    oscillatordataviewer.h
    oscillatordataviewer.cpp
)

target_include_directories(DrawSound PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DrawSound
    DiagWindow
<<<<<<< HEAD
    Qt6::Widgets
    Qt6::Multimedia
=======
    Qt5::Widgets
    Qt5::Multimedia
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
