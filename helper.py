from flask import Blueprint
from pymongo import MongoClient # Database connector
import pickle
import numpy as np

# This file contains helper methods

helper = Blueprint('helper', __name__)

client = MongoClient('localhost', 27017)    #Configure the connection to the database
db = client.data    #Select the database
dat = db.dictdata   #Select the collection

# loads and saves states and variables from the database into the 'states' dictionary
# TODO: add code to handle multiple requests
def load_states():
    if dat.find({ "_id": { "$exists": True, "$ne": "" } }):
        idinfo = dat.find_one({"_id":1})
        if 'states' in idinfo:
            #print("db", file=sys.stderr)
            states = pickle.loads(idinfo['states'])
            return states
        else:
            dat.update_one({'_id':1},{'$set':{'states': {} }})

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
