COMPILER = /usr/local/opt/llvm/bin/clang++
#COMPILER = g++ -g -Wall -std=c++11 -stdlib=libc++
# Python and BoostPython links.
BOOSTHEADERS = -I/usr/local/Cellar/boost/1.70.0/include/boost
BOOSTLIBRARIES = -L/usr/local/Cellar/boost-python3/1.70.0_1/lib/
# This is the boost library we want to use, there are also libraries for multithreading etc. 
# All we do is find the file libboost_python3.a and link it by taking away the 'lib' and '.a'.
BOOSTLIB = -lboost_python37
PYTHONHEADERS = -I/Users/gurgentus/anaconda3/include/python3.7m
PYTHONLIBRARIES = -L/Users/gurgentus/anaconda3/lib
# Link to python3.5 library, the same as we did for boost.
PYTHONLIB = -lpython3.7m

# Eigen
EIGEN_LIB = -I/Users/gurgentus/Code/eigen_dir

# Collect links.
LIBRARIES = $(BOOSTLIBRARIES) $(PYTHONLIBRARIES) $(PYTHONLIB) $(BOOSTLIB)
HEADERS = $(BOOSTHEADERS) $(PYTHONHEADERS) $(EIGEN_LIB)

MYSOURCES= $(wildcard *.cpp)
MYOBJECTS= $(SOURCES:.cpp=.o)

# Build target.

TARGET = OrbitTransfer


# Node.o: Node.h
# 	$(COMPILER) -c Node.h

# FiniteDifferenceGrid.o: Node.o FiniteDifferenceGrid.cpp FiniteDifferenceGrid.h
# 	$(COMPILER) -c FiniteDifferenceGrid.cpp

# DifferentialSystem.o: FiniteDifferenceGrid.o DifferentialSystem.cpp DifferentialSystem.h
# 	$(COMPILER) -c DifferentialSystem.cpp

# AbstractDeSolver.o: FiniteDifferenceGrid.o DifferentialSystem.o AbstractDeSolver.cpp AbstractDeSolver.h
# 	$(COMPILER) -c AbstractDeSolver.cpp

# RungeKuttaSolver.o: RungeKuttaSolver.cpp RungeKuttaSolver.h
# 	$(COMPILER) -c RungeKuttaSolver.cpp

# CollocationSolver.o: AbstractDeSolver.o CollocationSolver.cpp CollocationSolver.h
# 	$(COMPILER) -c CollocationSolver.cpp

# BEGIN MAKE
all: $(TARGET) Car

# set(SOURCE_FILES DifferentialSystem.cpp DifferentialSystem.h FiniteDifferenceGrid.cpp FiniteDifferenceGrid.h Node.h
#         AbstractDeSolver.cpp AbstractDeSolver.h CollocationSolver.cpp CollocationSolver.h
#         RungeKuttaSolver.cpp RungeKuttaSolver.h Car.cpp Car.h main.cpp


$(TARGET): $(TARGET).cpp
#  Node.h FiniteDifferenceGrid.cpp FiniteDifferenceGrid.h DifferentialSystem.cpp DifferentialSystem.h AbstractDeSolver.cpp AbstractDeSolver.h RungeKuttaSolver.cpp RungeKuttaSolver.h CollocationSolver.cpp CollocationSolver.h 
#  $(TARGET).cpp 
	# Note that '-shared' creates a library that is accessible.
	$(COMPILER) -shared -fPIC $(LIBRARIES) $(HEADERS) $(TARGET).cpp FiniteDifferenceGrid.cpp AbstractDeSolver.cpp DifferentialSystem.cpp RungeKuttaSolver.cpp CollocationSolver.cpp -o $(TARGET).so

Car: Car.cpp
	# Note that '-shared' creates a library that is accessible.
	$(COMPILER) -shared -fPIC $(LIBRARIES) $(HEADERS) $(TARGET).cpp FiniteDifferenceGrid.cpp AbstractDeSolver.cpp DifferentialSystem.cpp RungeKuttaSolver.cpp CollocationSolver.cpp  Car.cpp -o Car.so


clean:
	$(RM) $(TARGET)


# # location of the Python header files

# PYTHON_VERSION = 3.7m
# #PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
# PYTHON_INCLUDE = /Users/gurgentus/anaconda3/include/python3.7m
# PYTHON_LIB = /Users/gurgentus/anaconda3/lib/

# # Eigen
# EIGEN_LIB = /Users/gurgentus/Code/eigen_dir


# # location of the Boost Python include files and library

# #BOOST_INC = /usr/include

# BOOST_INC = /usr/local/Cellar/boost-python3/1.70.0_1/include
# #BOOST_LIB = /usr/lib
# BOOST_LIB = /usr/local/Cellar/boost/1.70.0/lib
# BOOSTP_LIB = /usr/local/Cellar/boost-python3/1.70.0_1/lib
# # compile mesh classes
# TARGET = omt

# $(TARGET).a: $(TARGET).o
# 	#g++ -shared $(TARGET).o -L$(BOOSTP_LIB) -lboost_python3 -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) -o $(TARGET).so
# 	g++ -shared -Wl -o $(TARGET).so -L$(BOOST_LIB) -L$(BOOSTP_LIB) -I$(EIGEN_LIB) -L$(PYTHON_LIB) -lpthread
# 	#g++ -o $(TARGET) -L$(BOOSTP_LIB) -lboost_python3 -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) -lpthread

# $(TARGET).o: $(TARGET).cpp
# 	g++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -I$(EIGEN_LIB) -fPIC -c $(TARGET).cpp
# #g++ -fPIC -I/usr/include/python3.6m library.cpp example.cpp -shared -o example.so -lboost_python3


# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
