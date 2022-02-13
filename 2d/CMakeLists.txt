
cmake_minimum_required(VERSION 3.1)
project(WaterSim)

# Finds Libigl
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

# Custom C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -Wall -Wextra -Wpedantic")
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

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/include/*.cpp)
file(GLOB HFILES *.h)
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_definitions(-DIGL_VIEWER_VIEWER_QUIET)
add_executable(${PROJECT_NAME} ${SRCFILES} ${SRC_FILES} ${HFILES} ${H_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png)
