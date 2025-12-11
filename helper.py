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
import os
from dotenv import load_dotenv

# This file contains helper methods

helper = Blueprint('helper', __name__)

# Load environment variables from .env file
load_dotenv()

# MongoDB connection - check once at startup
_use_mongodb = False
_client = None
_db = None
_dat = None

def _init_mongodb():
    """Initialize MongoDB connection once at startup"""
    global _use_mongodb, _client, _db, _dat

    mongodb_uri = os.getenv('MONGODB_URI', 'mongodb://localhost:27017/')
    try:
        # Try to connect with a short timeout
        _client = MongoClient(mongodb_uri, server_api=ServerApi('1'),
                            serverSelectionTimeoutMS=1000,
                            connectTimeoutMS=1000)
        # Test connection
        _client.admin.command('ping')
        _db = _client.data
        _dat = _db.dictdata
        _use_mongodb = True
        print("✅ Connected to MongoDB")
    except Exception as e:
        print(f"⚠️  MongoDB not available: {e}")
        print("📝 Using in-memory state storage")
        _use_mongodb = False
        _client = None
        _db = None
        _dat = None

# Initialize MongoDB connection once at module import
_init_mongodb()

# Fallback in-memory storage when MongoDB is not available
_memory_storage = {'states': {}}

# loads and saves states and variables from the database into the 'states' dictionary
# TODO: add code to handle multiple requests
def load_states():
    """Load states from MongoDB or memory fallback"""
    # Use the startup flag to decide - no repeated connection attempts!
    if not _use_mongodb:
        return _memory_storage.get('states', {})

    try:
        if _dat.find({ "_id": { "$exists": True, "$ne": "" } }):
            idinfo = _dat.find_one({"_id":1})
        else:
            _dat.update_one({"_id":1}, {"$set": {"id":1}}, upsert=True)
            idinfo = _dat.find_one({"_id":1})

        if idinfo and 'states' in idinfo:
            states = pickle.loads(idinfo['states'])
        else:
            _dat.update_one({'_id':1},{'$set':{'states': {} }}, upsert=True)
            idinfo = _dat.find_one({"_id":1})
            states = pickle.loads(idinfo['states'])
        return states
    except Exception as e:
        print(f"Warning: Could not load from MongoDB: {e}. Using in-memory storage.")
        return _memory_storage.get('states', {})

def save_states(states):
    """Save states to MongoDB or memory fallback"""
    # Use the startup flag to decide - no repeated connection attempts!
    if not _use_mongodb:
        _memory_storage['states'] = states
        return

    try:
        _dat.update_one({'_id':1},{'$set':{'states': pickle.dumps(states)}})
    except Exception as e:
        print(f"Warning: Could not save to MongoDB: {e}. Using in-memory storage.")
        _memory_storage['states'] = states

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
                   x_axis_label=x_label, y_axis_label=y_label, width=600, height=height)

    # create a color iterator
    #colors = itertools.cycle(palette)
    num_colors = len(x)
    # add a line renderer with legend and line thickness
    p.multi_line(x, y, legend_label=legend, color=bokeh.palettes.viridis(num_colors), line_width=2)

    script, div = components(p)
    states = load_states()
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}
    save_states(states)

def create_plot(x, y, plot, plot_title, x_label, y_label, legend, height=600):
    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title=plot_title, tools=TOOLS,
                   x_axis_label=x_label, y_axis_label=y_label, width=600, height=height)


    #ts = [2,3]
    #ys = [5,6]
    # add a line renderer with legend and line thickness
    p.line(x, y, legend_label=legend, line_width=2)

    script, div = components(p)
    states = load_states()
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}
    save_states(states)
