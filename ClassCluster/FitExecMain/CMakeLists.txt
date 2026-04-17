<<<<<<< HEAD
cmake_minimum_required(VERSION 3.5.0)
=======
cmake_minimum_required(VERSION 3.1.0)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

project(FitExecMain)

add_definitions(-DDATA_DIRECTORY="${CMAKE_SOURCE_DIR}/ClassCluster/Data")

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_executable(FitExecMain
    datasortfunctor.h
    datasortfunctor.cpp
    fitexeccontrol.h
    fitexeccontrol.cpp
    main.cpp
)

target_link_libraries(FitExecMain
    ControlWindow
    FitObject
	Logging
    MolSpektAnalysisMainWindow
    ParticleWatch
    Simulation
    TestHelper
<<<<<<< HEAD
    Qt6::Widgets
=======
    Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
