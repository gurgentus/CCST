# import standard packages
import os, math, pickle, json, sys
import numpy as np
import bokeh.plotting as plt
from scipy import linalg as la
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

UPLOAD_FOLDER = os.path.basename('uploads')
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
    else:
        return str(states[name]['value'].tolist())

@api.dispatcher.add_method
def orb_transfer_lobatto3(mu, m0, Isp, T, r0, tf, plot):
    orb_ins = OrbitTransfer.OrbitTransfer()
    status = orb_ins.run(float(mu), float(m0), float(Isp), float(T), float(r0), float(tf))
    if (status == 1):
        return 'Max number of iterations exhausted.'
    
    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Orbit Transfer", tools=TOOLS,
                   x_axis_label='x [km]', y_axis_label='y [km]', plot_width=800, plot_height=800)
    ts = orb_ins.getX()
    ys = orb_ins.getY()
    
    #ts = [2,3]
    #ys = [5,6]
    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend='r', line_width=2)
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}
    helper.save_states(states)
    return 'Orbital transfer plot created. Type gdisplay ' + plot + ' to display'

@api.dispatcher.add_method
def orb_transfer(mu, m0, Isp, T, r0, tf, k, rho, a, b, plot):
    orb_ins = OrbitTransfer.OrbitTransfer()
    orb_ins.SetMatrix(int(k), json.loads(rho), json.loads(a), json.loads(b))
    status = orb_ins.run(float(mu), float(m0), float(Isp), float(T), float(r0), float(tf))
    if (status == 1):
        return 'Max number of iterations exhausted.'

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Orbit Transfer", tools=TOOLS,
                   x_axis_label='x [km]', y_axis_label='y [km]', plot_width=800, plot_height=800)

    ts = orb_ins.getX()
    ys = orb_ins.getY()
    #ts = [2,3]
    #ys = [5,6]
    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend='r', line_width=2)

    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states(states)

    return 'Orbital transfer plot created. Type gdisplay ' + plot + ' to display'

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
