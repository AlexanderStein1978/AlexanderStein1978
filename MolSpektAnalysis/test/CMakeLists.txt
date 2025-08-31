cmake_minimum_required(VERSION 3.1.0)

cmake_policy(SET CMP0057 NEW)

project(MolSpektAnalysis_test)

include(GoogleTest)

find_package(Qt5 COMPONENTS Core REQUIRED)

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
	Qt5::Core
)

gtest_add_tests(TARGET MolSpektAnalysis_test)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/TestSpline.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/TestSpline2.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
