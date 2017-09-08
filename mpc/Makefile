# location of the Python header files

PYTHON_VERSION = 3.5m
#PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
PYTHON_INCLUDE = /Applications/anaconda/include/python3.5m
PYTHON_LIB = /Applications/anaconda/lib/

# Eigen
EIGEN_LIB = /Users/gurgentus/Dropbox/Code/ccst/Eigen-3.3


# location of the Boost Python include files and library

#BOOST_INC = /usr/include
BOOST_INC = /usr/local/Cellar/boost/1.65.0/include
#BOOST_LIB = /usr/lib
BOOST_LIB = /usr/local/Cellar/boost/1.65.0/lib
BOOSTP_LIB = /usr/local/Cellar/boost-python/1.65.0/lib
# compile mesh classes
TARGET = omt

$(TARGET).so: $(TARGET).o
	g++ -shared $(TARGET).o -L$(BOOSTP_LIB) -lboost_python3 -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) -o $(TARGET).so

$(TARGET).o: $(TARGET).cpp
	g++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -I$(EIGEN_LIB) -fPIC -c $(TARGET).cpp
