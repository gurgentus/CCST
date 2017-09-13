Heroku Flask Setup and Deployment Example:

http://blog.shea.io/lightweight-python-apps-with-flask-twitter-bootstrap-and-heroku/

Installing virtualenv with conda installed:

conda install virtualenv
pip install numpy
pip install json-rpc

Mongo installation:

https://treehouse.github.io/installation-guides/mac/mongo-mac.html

Boost.NumPy

$ git clone https://github.com/ndarray/Boost.NumPy.git
$ cd Boost.NumPy && mkdir build && cd build

In CMakeLists.txt change Python3 to Python-py35 or 34

$ cmake -DPYTHON_LIBRARY=$HOME/anaconda3/lib/libpython3.5m.so ../
$ make install

or

cmake -D Boost_NO_BOOST_CMAKE=ON -DCMAKE_INSTALL_PREFIX:PATH=/home/gurgentus/ccst ..

to build to a specific directory
