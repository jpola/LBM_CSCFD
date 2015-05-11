cmake_minimum_required(VERSION 2.8)

project(LBM_CPU)
set(APP LBM_CPU)

#provide path to openframeworks CMAKE
#(content of of dir : addons dev lib-linux  openFrameworks.cmake  src CMakeLists.txt  examples  licenses)

list(APPEND CMAKE_MODULE_PATH /home/jpola/Projects/OpenFrameworks/of)

include(openFrameworks)
ofxaddon(ofxGui)

#find OpenMP for parallel operations
find_package( OpenMP REQUIRED )
message(STATUS "OpenMP_CXX_FLAGS: " ${OpenMP_CXX_FLAGS})
message(STATUS "OpenMP_C_FLAGS: " ${OpenMP_C_FLAGS})
message(STATUS "Release Flags" ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "Debug Flags" ${CMAKE_CXX_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})


include_directories(src/colormap)
#include_directories()

set(INCLUDES
    src/node.hpp
    src/lbmTypes.hpp
    src/lbmController.hpp
)

set(SOURCES
    src/main.cpp
    src/colormap/ofxColorMap.cpp
    src/colormap/namedMaps.cpp
    src/ofApp.cpp
)

add_executable(
    ${APP}
    ${SOURCES}
    ${INCLUDES}
)

target_link_libraries(
    ${APP}
    ${OFXADDONS_BEGIN}
    ${OFXADDONS_LIBRARIES}
    ${OFXADDONS_END}
    ${OPENFRAMEWORKS_LIBRARIES}
)

#cotire(${APP})
