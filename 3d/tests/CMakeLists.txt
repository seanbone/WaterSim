cmake_minimum_required(VERSION 3.1)
project(WaterSim-Tests)

add_definitions(-DCMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR})

# Add custom targets
add_custom_target(watersim-tests-build)
add_custom_target(watersim-tests COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
add_dependencies(watersim-tests watersim-tests-build)

# Auto-detect files in this directory
file(GLOB TEST_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Make a test for each file
foreach(file ${TEST_FILES})
    get_filename_component(filename ${file} NAME_WE)
    message(STATUS "Adding test " ${filename} "\t" ${file})

    add_executable(${filename} ${file} ${PROJECT_SOURCE_DIR}/../src/NcReader.cpp)
    target_include_directories(${filename} PRIVATE includes)
    target_link_libraries(${filename} watersim-core netcdf-cxx4)

    add_test(NAME ${filename} COMMAND ${filename})
    add_dependencies(watersim-tests-build ${filename})
endforeach()

enable_testing()