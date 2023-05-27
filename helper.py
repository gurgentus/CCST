from flask import Blueprint
from pymongo.mongo_client import MongoClient
from pymongo.server_api import ServerApi
import pickle
import numpy as np
import bokeh.plotting as plt
from bokeh.resources import CDN
from bokeh.embed import components
import bokeh.palettes
import itertools

# This file contains helper methods

helper = Blueprint('helper', __name__)

uri = 'mongodb+srv://gurgentus:e5O2rrDipb9nm5Wj@ccst.cf030e2.mongodb.net/?retryWrites=true&w=majority'
#uri = 'mongodb://heroku_c9chv2pq:euoqe7c7o24l17ar4pavqleame@.mlab.com:21014/heroku_c9chv2pq'

#client = MongoClient('localhost', 27017)    #Configure the connection to the database
# client = MongoClient(uri)
# Create a new client and connect to the server
client = MongoClient(uri, server_api=ServerApi('1'))

db = client.data    #Select the database
#db = client.get_default_database()

dat = db.dictdata   #Select the collection

# loads and saves states and variables from the database into the 'states' dictionary
# TODO: add code to handle multiple requests
def load_states():
    if dat.find({ "_id": { "$exists": True, "$ne": "" } }):
        idinfo = dat.find_one({"_id":1})
    else:
        dat.update_one({"_id":1}, {"$set": {"id":1}}, upsert=True)
        idinfo = dat.find_one({"_id":1})

    if idinfo and 'states' in idinfo:
        #print("db", file=sys.stderr)
        states = pickle.loads(idinfo['states'])
    else:
        dat.update_one({'_id':1},{'$set':{'states': {} }}, upsert=True)

    idinfo = dat.find_one({"_id":1})
    states = pickle.loads(idinfo['states'])
    return states

def save_states(states):
    dat.update_one({'_id':1},{'$set':{'states': pickle.dumps(states)}})

# represent the matrix as a string
def to_str_repr(arg):
    x = np.array(arg)
    y = np.array([str(w) for w in x.reshape(x.size)])
    return y.reshape(arg.shape).tolist()

# represent a numpy array as string
def nparray_to_str_repr(x):
    y = np.array([str(w) for w in x.reshape(x.size)])
    return y.reshape(x.shape).tolist()

def create_multi_plot(x, y, plot, plot_title, x_label, y_label, legend, height=600):
    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title=plot_title, tools=TOOLS,
                   x_axis_label=x_label, y_axis_label=y_label, plot_width=600, plot_height=height)

    # create a color iterator
    #colors = itertools.cycle(palette)
    num_colors = len(x)
    # add a line renderer with legend and line thickness
    p.multi_line(x, y, legend=legend, color=bokeh.palettes.viridis(num_colors), line_width=2)

    script, div = components(p)
    states = load_states()
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}
    save_states(states)

def create_plot(x, y, plot, plot_title, x_label, y_label, legend, height=600):
    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title=plot_title, tools=TOOLS,
                   x_axis_label=x_label, y_axis_label=y_label, plot_width=600, plot_height=height)


    #ts = [2,3]
    #ys = [5,6]
    # add a line renderer with legend and line thickness
    p.line(x, y, legend=legend, line_width=2)

    script, div = components(p)
    states = load_states()
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}
    save_states(states)
