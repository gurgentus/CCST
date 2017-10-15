cmake_minimum_required(VERSION 2.8)
project(OrbitTransfer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_STANDARD 11)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)
##FIND_PACKAGE(Boost)
FIND_PACKAGE(Boost COMPONENTS python-py34)
#
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} "/home/gurgentus/ccst/include")
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

set(SOURCE_FILES DifferentialSystem.cpp DifferentialSystem.h FiniteDifferenceGrid.cpp FiniteDifferenceGrid.h Node.h
        AbstractDeSolver.cpp AbstractDeSolver.h CollocationSolver.cpp CollocationSolver.h
        RungeKuttaSolver.cpp RungeKuttaSolver.h Car.cpp Car.h main.cpp)
#add_executable(OrbitalTransfer ${SOURCE_FILES})
add_library(OrbitTransfer SHARED ${SOURCE_FILES})
target_link_libraries(OrbitTransfer ${Python_LIBRARIES} ${Boost_LIBRARIES})
add_library(Car SHARED ${SOURCE_FILES})
target_link_libraries(Car ${Python_LIBRARIES} ${Boost_LIBRARIES})

