cmake_minimum_required(VERSION 3.8)
project(nums)

set(CMAKE_CXX_STANDARD 11)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(Boost COMPONENTS python3)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

set(SOURCE_FILES main.cpp BVPSolver.cpp BVPSolver.h DifferentialSystem.cpp DifferentialSystem.h FiniteDifferenceGrid.cpp FiniteDifferenceGrid.h Node.h OrbitTransfer.cpp OrbitTransfer.h CollocationSolver.cpp CollocationSolver.h)# PythonWrapper.cpp)
#add_executable(nums ${SOURCE_FILES})
ADD_LIBRARY(nums SHARED ${SOURCE_FILES})
target_link_libraries(nums ${Python_LIBRARIES} ${Boost_LIBRARIES})


