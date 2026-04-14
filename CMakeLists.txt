# CMakeLists.txt
# cmake -B build -DCMAKE_BUILD_TYPE=Debug
# cmake --build build --parallel 4
# cmake --build build -j 8
# cd build: ctest --output-on-failure

cmake_minimum_required(VERSION 3.15)

if(DEFINED ENV{CI_VERSION})
    project(RelSpacetime VERSION $ENV{CI_VERSION})
    message(STATUS "Building CI Version: ${PROJECT_VERSION}")
else()
    project(RelSpacetime VERSION 1.0.0) # FALLBACK
    message(STATUS "Building Local Version: ${PROJECT_VERSION}")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0 
)

# Download Eigen but do not evaluate its CMakeLists.txt
FetchContent_GetProperties(Eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(Eigen)
    
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE ${eigen_SOURCE_DIR})
endif()

add_library(RelSpacetimeCore
    src/AnalyticalRaytracer.cpp
    src/DelayedTimeSolver.cpp
    src/MetricTensor.cpp
    src/SpacetimeVector.cpp
    src/LocalFrame.cpp
    src/WorldlineIntegrator.cpp
    src/Transform.cpp
    src/SpacetimeInterpolator.cpp
    src/VisualFX.cpp
)

target_include_directories(RelSpacetimeCore PUBLIC include)

target_link_libraries(RelSpacetimeCore PUBLIC Eigen3::Eigen)

set_target_properties(RelSpacetimeCore PROPERTIES POSITION_INDEPENDENT_CODE ON)

# CTest
enable_testing()
add_subdirectory(tests)

# GODOT GDEXTENSION SETUP

# 1. Tell godot-cpp where your custom API file is BEFORE adding the subdirectory
set(GODOT_CUSTOM_API_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/extension_api.json" CACHE STRING "Path to custom extension_api.json" FORCE)

# 2. Add godot-cpp. CMake will automatically look for the CMakeLists.txt inside this folder!
add_subdirectory(ThirdParty/godot-cpp)

# 3. Create your Godot extension shared library (.dylib)
add_library(RelSpacetimeGodot SHARED
    src_godot/register_types.cpp
    src_godot/GDEntityState.cpp
    src_godot/GDRelativityAPI.cpp
)

# 4. Link your wrapper to BOTH your core physics SDK and godot-cpp
target_link_libraries(RelSpacetimeGodot 
    PRIVATE 
        RelSpacetimeCore 
        godot-cpp
)

# 5. Output the .dylib directly into your bin/ folder so Godot can find it
set_target_properties(RelSpacetimeGodot PROPERTIES
    PREFIX "lib"
    SUFFIX ".macos.debug.framework.arm64.dylib" 
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../Documents/rel-spacetime-test/bin"
)

