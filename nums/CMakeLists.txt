cmake_minimum_required(VERSION 3.11.0)
 
project(OrbitTransfer)
# Find necessary packages
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIR})
 
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#add_definitions     ( ${EIGEN3_DEFINITIONS} )
#include_directories ( ${EIGEN3_INCLUDE_DIRS} )
#include_directories(${CMAKE_SOURCE_DIR}/)
include_directories(/Users/gurgentus/Code/CCST/nums/CMakeLists.txt)
#find_package(Boost COMPONENTS python3 REQUIRED)

find_package(Boost COMPONENTS system filesystem REQUIRED)

#FIND_PACKAGE(Boost COMPONENTS python-py37m)

include_directories(${Boost_INCLUDE_DIR})
 
# Build & Link our library
#add_library(hello MODULE hello.cpp)

if(APPLE)
        set(PYTHON_LIBRARIES "/usr/local/Cellar/python/3.7.3/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib")
endif()

set (CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES DifferentialSystem.cpp DifferentialSystem.h FiniteDifferenceGrid.cpp FiniteDifferenceGrid.h Node.h
        AbstractDeSolver.cpp AbstractDeSolver.h CollocationSolver.cpp CollocationSolver.h
        RungeKuttaSolver.cpp RungeKuttaSolver.h Car.cpp Car.h main.cpp)
#add_executable(OrbitalTransfer ${SOURCE_FILES})
add_library(OrbitTransfer SHARED ${SOURCE_FILES})
target_link_libraries(OrbitTransfer ${PYTHON_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
add_library(Car SHARED ${SOURCE_FILES})
target_link_libraries(Car ${PYTHON_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# don't prepend wrapper library name with lib
#set_target_properties(OrbitTransfer PROPERTIES PREFIX "")