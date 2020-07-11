cmake_minimum_required(VERSION 3.1.0)

include(GoogleTest)

project(MolSpektAnalysis_test)

add_executable(MolSpektAnalysis_test
	integrationtest.h
	integrationtest.cpp
	test_potential.cpp
	main.cpp
)

target_link_libraries(MolSpektAnalysis_test
	gtest
	Potential
)

gtest_add_tests(TARGET MolSpektAnalysis_test)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/TestSpline.pot DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
