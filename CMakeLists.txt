cmake_minimum_required(VERSION 3.1.0)
project(Physics)

cmake_policy(SET CMP0057 NEW)

set(CMAKE_CXX_FLAGS_RELEASE "-Ofast" CACHE STRING "" FORCE)

add_subdirectory(ClassCluster)
add_subdirectory(MolSpektAnalysis)
add_subdirectory(Sound)
add_subdirectory(googletest)
