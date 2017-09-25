# import standard packages
import os, math, pickle, json, sys
import numpy as np
import bokeh.plotting as plt
from scipy import linalg as la
import scipy.misc
from PIL import Image
from scipy.integrate import ode
from flask import Flask, jsonify, render_template, send_from_directory, flash, request, redirect, url_for
from jsonrpc.backend.flask import api
from pymongo import MongoClient # Database connector
from bson.objectid import ObjectId # For ObjectId to work
from bson.binary import Binary
# from oct2py import octave as oc
from data import data_api
from control import control_api
from flask_jsglue import JSGlue

from bokeh.resources import CDN
from bokeh.embed import components

# orbital mechanics toolbox
import mpc.omt as omt
import nums.OrbitTransfer as OrbitTransfer
# loads data from database and contains helper methods
import helper


# flask initialization
app = Flask(__name__)
app.config.update(
    DEBUG = True,
)

jsglue = JSGlue(app)

#app.add_url_rule('/', 'api', api.as_view())
app.register_blueprint(api.as_blueprint())
app.register_blueprint(control_api)
app.register_blueprint(data_api)

#UPLOAD_FOLDER = os.path.basename('uploads')
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'uploads')

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route("/")
def index():
    return render_template('calcresult.html', states=states)

# controllers

@app.route('/upload', methods=['POST'])
def upload_file():
    file = request.files['file']
    f = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
    file.save(f)
    name = request.form['inputTitle']
    str_value = file.filename

    states = helper.load_states()
    states[name] = {'value': str_value, 'meta': {'what': 'image', 'value': str_value}}
    helper.save_states(states)
    flash("File Uploaded.  To show the image type: gdisplay " + name)
    return render_template('calcresult.html', img_name=file.filename)

@app.route('/uploads/<filename>')
def send_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)

@api.dispatcher.add_method
def generate_grid(name, height, width):
    image_shape = (height, width)
    states = helper.load_states()
    f = os.path.join(app.config['UPLOAD_FOLDER'], states[name]['value'])
    image = scipy.misc.imresize(scipy.misc.imread(name=f, flatten=True), image_shape)
    str_vec = helper.to_str_repr(np.array(image))
    states[name+'_grid'] = {'value': np.matrix(image), 'meta': {'what': 'matrix', 'value': str_vec}}
    helper.save_states(states)
    #matrix(name+'_grid', np.array(image))
    return 'Grid matrix generated.  Type gdisplay ' + name + '_grid to display.'

@api.dispatcher.add_method
def generate_grid_obstacles(name, threshold_value, new_name):
    states = helper.load_states()
    array_np = np.array(states[name]['value'])
    more_then = array_np >= threshold_value
    less_then = array_np < threshold_value
    array_np[more_then] = 1
    array_np[less_then] = 0
    # print(array_np, file=sys.stderr)
    str_vec = helper.to_str_repr(array_np)
    states[new_name] = {'value': np.matrix(array_np), 'meta': {'what': 'matrix', 'value': str_vec}}
    helper.save_states(states)
    return 'New grid matrix generated.  Type gdisplay ' + new_name + ' to display.'

@api.dispatcher.add_method
def dijkstras(name, start, goal, result): #,x_spacing,y_spacing,start,goal):

    # Implements Dijkstra's shortest path algorithm
    # Input:
    # occupancyshhaa_map - an N by M numpy array of boolean values (represented
    #     as integers 0 and 1) that represents the locations of the obstacles
    #     in the world
    # x_spacing - parameter representing spacing between adjacent columns
    # y_spacing - parameter representing spacing between adjacent rows
    # start - a 3 by 1 numpy array of (x,y,theta) for the starting position
    # goal - a 3 by 1 numpy array of (x,y,theta) for the finishing position
    # Output:
    # path: list of the indices of the nodes on the shortest path found
    #     starting with "start" and ending with "end" (each node is in
    #     metric coordinates)

    states = helper.load_states()
    occupancy_map = np.array(states[name]['value'])
    color_map = np.array(states[name]['value'])
    # color_map = np.asmatrix(color_map)
    # print(color_map, file=sys.stderr)
    # print(occupancy_map, file=sys.stderr)

    INF = math.inf #100000
    mymap = np.zeros_like(occupancy_map)

    nrows = np.shape(mymap)[0]
    ncols = np.shape(mymap)[1]
    #print(ncols,nrows)
    start = np.array(json.loads(start))
    goal = np.array(json.loads(goal))

    start_indx = start[1] #int(np.ceil((start[0]/x_spacing)-0.5))
    start_indy = start[0] #int(np.ceil((start[1]/y_spacing)-0.5))
    dest_indx = goal[1] #int(np.ceil((goal[0]/x_spacing)-0.5))
    dest_indy = goal[0] #int(np.ceil((goal[1]/y_spacing)-0.5))

    distanceFromStart = np.zeros_like(occupancy_map, float)
    distanceFromStart.fill(INF)
    parent = np.zeros_like(occupancy_map)

    mymap[np.where(occupancy_map == 0)] = 1  # Mark free cells
    mymap[np.where(occupancy_map == 1)] = 2   # Mark obstacle cells

    # Generate linear indices of start and dest nodes
    dest_node = np.ravel_multi_index([dest_indy, dest_indx], np.shape(mymap), order='C')

    mymap[start_indy][start_indx] = 5
    mymap[dest_indy][dest_indx]  = 6

    distanceFromStart[start_indy][start_indx] = 0

    # keep track of number of nodes expanded
    numExpanded = 0
    # Main Loop
    count = 0
    while True:
        count = count + 1
        if (count == 1000):
            return "Iteration limit exceeded."

        # Find the node with the minimum distance
        min_dist = np.min(distanceFromStart)
        current = np.argmin(distanceFromStart)
        # Compute row, column coordinates of current node
        [i,j] = np.unravel_index(current, np.shape(mymap), order='C')

        if (current == dest_node or min_dist == INF):
            break

        # Update map
        mymap[i][j] = 3         # mark current node as visited

        distanceFromStart[i][j] = INF; # remove this node from further consideration
        numExpanded = numExpanded + 1;
        # Visit each neighbor of the current node and update the map, distances
        # and parent tables appropriately.

        if (i>0) and (mymap[i-1][j] != 3) and (mymap[i-1][j] != 5) and (mymap[i-1][j] != 2):
            if (min_dist+1 < distanceFromStart[i-1][j]):
                distanceFromStart[i-1][j] = min_dist+1
                mymap[i-1][j] = 4
                parent[i-1][j] = current

        if (i<nrows-1) and (mymap[i+1][j] != 3) and (mymap[i+1][j] != 5) and (mymap[i+1][j] != 2):
            if (min_dist+1 < distanceFromStart[i+1][j]):
                distanceFromStart[i+1][j] = min_dist+1;
                mymap[i+1][j]=4
                parent[i+1][j] = current

        if (j>0) and (mymap[i][j-1] != 3) and (mymap[i][j-1] != 5) and (mymap[i][j-1] != 2):
            if (min_dist+1 < distanceFromStart[i][j-1]):
                distanceFromStart[i][j-1] = min_dist+1
                mymap[i][j-1]=4
                parent[i][j-1] = current

        if (j<ncols-1) and (mymap[i][j+1] != 3) and (mymap[i][j+1] != 5) and (mymap[i][j+1] != 2):
            if (min_dist+1 < distanceFromStart[i][j+1]):
                distanceFromStart[i][j+1] = min_dist+1;
                mymap[i][j+1]=4
                parent[i][j+1] = current

        #print(distanceFromStart, file=sys.stderr)
        #print(mymap, file=sys.stderr)

    # Construct route from start to dest by following the parent links
    if (distanceFromStart[dest_indy][dest_indx] == INF):
        routep = np.array([])
    else:
        routep = np.array([[dest_indy,dest_indx]])

        [i,j] = np.unravel_index(dest_node, np.shape(mymap))
        while (parent[i][j] != 0):
            [i,j] = np.unravel_index(parent[i][j], np.shape(mymap))
            routep = np.concatenate((routep, np.array([[i,j]])), axis=0)
            #print(i,j)
        #[i,j] = np.unravel_index(start_node, np.shape(mymap))
        #routep = np.concatenate((routep, np.array([[i,j]])), axis=0)

    s = len(routep)
    a = start[0]#[0]
    b = start[1]#[0]

    route = np.array([[a,b]])

    for k in range(1,s+1):

        [a,b]= (routep[s-k])#+0.5)*[x_spacing, y_spacing]
        color_map[a][b] = 2
        route = np.concatenate((route, np.array([[a,b]])), axis=0)

    a = goal[0]#[0]
    b = goal[1]#[0]
    route = np.concatenate((route, np.array([[a,b]])), axis=0)

    str_vec = helper.to_str_repr(color_map)

    states[result] = {'value': color_map, 'meta': {'what': 'color_matrix', 'value': str_vec}}
    helper.save_states(states)

    return route.tolist()
    #pass

# simple test method - ping pong
@api.dispatcher.add_method
def ping():
    return 'pong'

@api.dispatcher.add_method
def start():
    dat.insert_one({
        '_id':1 ,'states':states
        })

states = helper.load_states()

@api.dispatcher.add_method
def display(name):
    #print("test", file=sys.stderr)
    states = helper.load_states()
    if states[name]['meta']['what'] == 'scalar':
        return states[name]['value']
    elif states[name]['meta']['what'] == 'image':
        return "To show the image type: " + states[name]['meta']['value']
    elif states[name]['meta']['what'] == 'matrix':
        return str(states[name]['value'].tolist())
    else:
        return str(states[name]['value'].tolist())

@api.dispatcher.add_method
def orb_transfer_lobatto3(mu, m0, Isp, T, r0, days, timestep_hrs, N, plot):
    orb_ins = OrbitTransfer.OrbitTransfer()
    status = orb_ins.run(float(mu), float(m0), float(Isp), float(T), float(r0), float(days), float(timestep_hrs), int(N))
    if (status == 1):
        return 'Max number of iterations exhausted.'

    if (status == 2):
        return 'Collocation scheme failed.  Try decreasing the timestep.'

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Orbit Transfer", tools=TOOLS,
                   x_axis_label='x [km]', y_axis_label='y [km]', plot_width=600, plot_height=600)
    ts = orb_ins.getX()
    ys = orb_ins.getY()
    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend='r', line_width=2)
    #p.circle(ts, ys, size=2, color="navy", alpha=0.5)
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Steering Control", tools=TOOLS,
                   x_axis_label='t [sec]', y_axis_label='angle [deg]', plot_width=600, plot_height=300)
    ts = orb_ins.getT()
    ys = orb_ins.getAngle()
    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend='steering', line_width=2)
    script, div = components(p)
    plot_name = plot + '_control'
    states[plot_name] = {'meta': {}}
    states[plot_name]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states(states)

    return 'Orbital transfer plot created. Type gdisplay ' + plot + ' to display the trajectory or gdisplay ' + plot + '_control to display the steering angle'


@api.dispatcher.add_method
def orb_transfer(mu, m0, Isp, T, r0, tf, timestep_hrs, N, k, rho, a, b, plot):
    orb_ins = OrbitTransfer.OrbitTransfer()
    orb_ins.SetMatrix(int(k), json.loads(rho), json.loads(a), json.loads(b))
    status = orb_ins.run(float(mu), float(m0), float(Isp), float(T), float(r0), float(tf), float(timestep_hrs), int(N))
    if (status == 1):
        return 'Max number of iterations exhausted.'

    if (status == 2):
        return 'Collocation scheme failed.  Try decreasing the timestep.'

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Orbit Transfer", tools=TOOLS,
                   x_axis_label='x [km]', y_axis_label='y [km]', plot_width=600, plot_height=600)

    ts = orb_ins.getX()
    ys = orb_ins.getY()
    #ts = [2,3]
    #ys = [5,6]
    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend='r', line_width=2)

    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Steering Control", tools=TOOLS,
                   x_axis_label='t [sec]', y_axis_label='angle [deg]', plot_width=600, plot_height=300)
    ts = orb_ins.getT()
    ys = orb_ins.getAngle()
    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend='steering', line_width=2)
    script, div = components(p)
    plot_name = plot + '_control'
    states[plot_name] = {'meta': {}}
    states[plot_name]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states(states)

    return 'Orbital transfer plot created. Type gdisplay ' + plot + ' to display the trajectory or gdisplay ' + plot + '_control to display the steering angle'


@api.dispatcher.add_method
def lambert(r1, r2, t, prograde, mu, name):
    omt_ins = omt.omt()
    omt_ins.lambert(json.loads(r1), json.loads(r2), float(t), bool(prograde), float(mu))
    #omt_ins.lambert([5000,10000,2100], [-14600,2500,7000], 3600, True, 398600)
    #return "Initial position: " + str(r1) + ", initial velocity: " + str(omt_ins.get_v0())
    states = helper.load_states()

    states[name] = {'value': [r1, omt_ins.get_v0()], 'meta': {'what': 'orbit', 'h': omt_ins.get_h(), 'a': omt_ins.get_a(),
      'e': omt_ins.get_e(), 'Omega': omt_ins.get_Omega(), 'i': omt_ins.get_i(), 'omega': omt_ins.get_omega()}}
    helper.save_states(states)

    return "Orbit calculated.  To display orbital elements type gdisplay " + name + "."

# retrieves the variable passed as the argument to render using latex (see calcresult.html)
@app.route('/gdisplay')
def gdisplay():
    states = helper.load_states()
    name = request.args.get('name')
    return jsonify(states[name]['meta'])

# parse the command as an octave command
# @app.route('/octave')
# def octave():
#     #states = load_states()
#     name = request.args.get('name')
#     sp = name.split('=')
#     # print(len(sp), file=sys.stderr)
#     if len(sp) > 1:
#         rhs = sp[1].strip()
#         nm = sp[0].strip()
#     else:
#         rhs = sp[0].strip()
#         nm = 'ans'
#
#     val = oc.eval(rhs)
#
#     if type(val) is float:
#         states[nm] = {'value': val, 'meta': {'what': 'scalar', 'value': str(val)}}
#     if type(val) is np.ndarray:
#         val = np.asmatrix(val)
#         str_vec = helper.to_str_repr(val)
#         states[nm] = {'value': val, 'meta': {'what': 'matrix', 'value': str_vec}}
#
#     helper.save_states(states)
#
#     return jsonify(states[nm]['meta'])


# @api.dispatcher.add_method
# def add_matrices(mat1, mat2):
#     return (states[mat1]+states[mat2]).tolist()
#
# @app.route('/_process')
# def _process():
#     a = request.args.get('q')
#     b_new = json.loads(a)
#     c_new = np.matrix(b_new)
#     result = json.dumps(c_new.tolist())
#     return jsonify(result2=result)
#
# @app.route('/add/<int:num1>/<int:num2>/')
# def add(num1, num2):
#     num = num1 + num2
#     return render_template('calcresult.html', result=num)
#
# @app.route("/sin/<int:num1>/")
# def sin(num1):
#     result = math.sin(num1)
#     result = np.matrix( num1 )
#     return render_template('calcresult.html', result=result)
#
# @app.route('/favicon.ico')
# def favicon():
#     return send_from_directory(os.path.join(app.root_path, 'static'), 'ico/favicon.ico')
#
# @app.route('/_add_numbers')
# def add_numbers():
#     a = request.args.get('a', 0, type=int)
#     b = request.args.get('b', 0, type=int)
#     return jsonify(result=a + b)

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

# launch
app.secret_key = os.urandom(24)
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
