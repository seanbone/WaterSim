
cmake_minimum_required(VERSION 3.1)
project(WaterSim)

option(FETCHCONTENT_QUIET OFF)

# Fetch JSON library
include(FetchContent)
FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
        GIT_TAG v3.9.1)
FetchContent_GetProperties(json)
if(NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Fetch NetCDF library
include(FetchContent)
FetchContent_Declare(
        netcdf
        GIT_REPOSITORY https://github.com/Unidata/netcdf-cxx4
        GIT_TAG 162d78a8c4dc33379df746dfe232d0aaebc705db)
FetchContent_GetProperties(netcdf)
if(NOT netcdf_POPULATED)
    FetchContent_Populate(netcdf)
    set(NCXX_ENABLE_TESTS OFF CACHE BOOL "Enable tests. Run with 'make test'." FORCE)
    add_subdirectory(${netcdf_SOURCE_DIR} ${netcdf_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()




# Finds Libigl
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

# Custom C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -Wall -Wextra -Wpedantic")
# Ignore some warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder -Wno-unknown-pragmas -Wno-unused-parameter -Wno-misleading-indentation")

# Export compilation commands for VIM editor
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

if (NOT LIBIGL_FOUND)
    find_package(LIBIGL REQUIRED QUIET)
endif()
# Add project files
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB SRC_FILES_CLI src/FLIP.cpp src/Mac3d.cpp src/MeshExporter.cpp src/Particle.cpp
         src/SimConfig.cpp src/WaterSim.cpp)
file(GLOB H_FILES_CLI src/FLIP.h src/Mac3d.h src/MeshExporter.h src/Particle.h
        src/SimConfig.h src/WaterSim.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(../lib/eigen3)


## Enable tests
add_subdirectory(tests)
enable_testing()


add_definitions(-DIGL_VIEWER_VIEWER_QUIET)
add_executable(watersim-gui watersim-gui.cpp ${SRC_FILES} ${H_FILES})
add_executable(watersim-cli watersim-cli.cpp ${SRC_FILES_CLI} ${H_FILES_CLI})
add_executable(viewmesh viewmesh.cpp)
add_executable(netcdf-test netcdf-test.cpp)
target_link_libraries(watersim-gui igl::core igl::opengl_glfw igl::opengl_glfw_imgui
        igl::png nlohmann_json::nlohmann_json netcdf-cxx4)
target_link_libraries(watersim-cli igl::core nlohmann_json::nlohmann_json netcdf-cxx4)
target_link_libraries(viewmesh igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png)
target_link_libraries(netcdf-test netcdf-cxx4)
