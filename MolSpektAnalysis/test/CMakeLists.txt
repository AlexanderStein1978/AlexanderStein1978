cmake_minimum_required(VERSION 3.5.0)

cmake_policy(SET CMP0057 NEW)

project(MolSpektAnalysis_test)

include(GoogleTest)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Core REQUIRED)
=======
find_package(Qt5 COMPONENTS Core REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_executable(MolSpektAnalysis_test
	integrationtest.h
	integrationtest.cpp
	test_potential.cpp
	main.cpp
)

target_link_libraries(MolSpektAnalysis_test
	gtest
	Calc
	LineTable
	MdiChild
	MonteCarloSim
	Potential
	PotentialStructs
	Residualfit
	TableWindow
	Tools
<<<<<<< HEAD
	Qt6::Core
=======
	Qt5::Core
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)

gtest_add_tests(TARGET MolSpektAnalysis_test)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/TestSpline.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/TestSpline2.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
