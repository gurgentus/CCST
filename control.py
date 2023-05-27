import numpy as np
import bokeh.plotting as plt
import csv, json, sys
from flask import Blueprint
from jsonrpc.backend.flask import api
from scipy import linalg as la
from scipy.integrate import ode
from bokeh.resources import CDN
from bokeh.embed import components
import helper

control_api = Blueprint('control_api', __name__)
# @simple_page.route('/<page>')
# def show(page):

# define a matrix
@api.dispatcher.add_method
def matrix(name, value):
    states = helper.load_states()
    val = np.matrix(json.loads(value))
    str_vec = helper.to_str_repr(val)
    #print(states, file=sys.stderr)

    states[name] = {'value': val, 'meta': {'what': 'matrix', 'value': str_vec}}
    helper.save_states(states)

    return states[name]['value'].tolist()

# calculate eigenvectors/eigenvalues of the matrix mat.
@api.dispatcher.add_method
def modes(mat, V,E):
    states = helper.load_states()

    e_vals, e_vecs = np.linalg.eig(states[mat]['value'])
    str_vec = helper.to_str_repr(e_vecs)
    str_vals = helper.to_str_repr(e_vals)
    states[V] = {'value': np.matrix(e_vecs), 'meta': {'what': 'matrix', 'value': str_vec}}
    states[E] = {'value': e_vals, 'meta': {'what': 'vector', 'value': str_vals}}

    helper.save_states(states)

    return states[E]['meta']['value']

# define state-space system using matrices A,B,C,D and store it under variable name
@api.dispatcher.add_method
def ss(A, B, C, D, name):
    states = helper.load_states()
    states[name] = {'value': [A, B, C, D], 'A': states[A], 'B': states[B], 'C': states[C],
     'D': states[D], 'controls': {}, 'outputs': {}, 'meta': {'what': 'ss',
     'A': states[A]['meta'], 'B': states[B]['meta'], 'C': states[C]['meta'], 'D': states[D]['meta']}}
    helper.save_states(states)
    result = 'x\' = ' + A + 'x + ' + B + 'u , y = ' + C + 'x + ' + D + 'u'
    return result

# define a controller using matrices A, B1, B2, C, D1, D2 and store under name
@api.dispatcher.add_method
def controller(A, B1, B2, C, D1, D2, name):
    states = helper.load_states()
    states[name] = {'value': [A, B1, B2, C, D1, D2], 'A': states[A], 'B1': states[B1], 'B2': states[B2],
     'C': states[C], 'D1': states[D1], 'D2': states[D2], 'controls': {}, 'outputs': {},
     'meta': {'what': 'controller',
     'A': states[A]['meta'], 'B1': states[B1]['meta'], 'B2': states[B2]['meta'], 'C': states[C]['meta'],
     'D1': states[D1]['meta'], 'D2': states[D2]['meta']}}
    helper.save_states(states)
    result = 'x_c\' = ' + A + 'x_c + ' + B1 + 'y + ' + B2 + 'r , u = ' + C + 'x_c + ' + D1 + 'y + ' + D2 + 'r'
    return result

# connect state space system G and a controller K into a feedback loop and store
# the resulting system under name
@api.dispatcher.add_method
def feedback(G, K, name):
    states = helper.load_states()

    Ap = states[G]['A']['value']
    Bp = states[G]['B']['value']
    Cp = states[G]['C']['value']
    Dp = states[G]['D']['value']
    Ac = states[K]['A']['value']
    Bc1 = states[K]['B1']['value']
    Bc2 = states[K]['B2']['value']
    Cc = states[K]['C']['value']
    Dc1 = states[K]['D1']['value']
    Dc2 = states[K]['D2']['value']

    I = np.identity(Dc1.shape[0])
    Z = I-Dc1*Dp
    Zinv = np.linalg.inv(Z)

    A11 = Ap + Bp*Zinv*Dc1*Cp
    A12 = Bp*Zinv*Cc
    I = np.identity(Dp.shape[0])
    A21 = Bc1*(I+Dp*Zinv*Dc1)*Cp
    A22 = Ac + Bc1*Dp*Zinv*Cc
    A = np.bmat([[A11, A12], [A21, A22]])
    A_newname = states[G]['value'][0]+'_cl'

    B11 = Bp*Zinv*Dc2
    B21 = Bc2 + Bc1*Dp*Zinv*Dc2
    B = np.bmat([[B11], [B21]])
    B_newname = states[G]['value'][1]+'_cl'

    C11 = I+Dp*Zinv*Dc1
    C12 = Cp*Dp*Zinv*Cc
    C = np.bmat([[C11, C12]])
    C_newname = states[G]['value'][2]+'_cl'

    D = Dp*Zinv*Dc2
    D_newname = states[G]['value'][3]+'_cl'

    helper.matrix(A_newname, A)
    helper.matrix(B_newname, B)
    helper.matrix(C_newname, C)
    helper.matrix(D_newname, D)

    states = helper.load_states()

    states[name] = {'value': [A_newname, B_newname, C_newname, D_newname],
     'A': states[A_newname], 'B': states[B_newname], 'C': states[C_newname], 'D': states[D_newname],
     'controls': {}, 'outputs': {}, 'meta': {'what': 'feedback', 'A': states[A_newname]['meta'],
     'B': states[B_newname]['meta'], 'C': states[C_newname]['meta'], 'D': states[D_newname]['meta']}}
    helper.save_states(states)
    result = 'x_a\' = ' + A_newname + 'x_a + ' + B_newname + 'r , y = ' + C_newname + 'x_a + ' + D_newname + 'r'
    return result

# associate the num's control variable with name
@api.dispatcher.add_method
def control(num, G, name):
    states = helper.load_states()
    states[G]['controls'][name] = num - 1
    helper.save_states(states)
    return 'Control number ' + str(num) + ' of ' + G + ' is set to ' + name

# associate the num's output variable with name
@api.dispatcher.add_method
def output(num, G, name):
    states = helper.load_states()
    states[G]['outputs'][name] = num - 1
    helper.save_states(states)
    return 'Output number ' + str(num) + ' of ' + G + ' is set to ' + name

# simulate output response of the state-space system G to initial conditions init
# and a constant forcing term inp for a time period time and save the
# resulting graph to plot
@api.dispatcher.add_method
def simulate(output, G, init, inp, time, plot):
    states = helper.load_states()

    t1 = time
    dt = 0.1
    ts = []
    ys = []

    def f(t, y, arg1):
        return states[G]['A']['value']*y.reshape(y.size,1) + states[G]['B']['value']*inp
    def jac(t, y, arg1):
        return states[G]['A']['value']

    def solout(t, y):
        ts.append(t)
        ys.append(y)

    r = ode(f, jac).set_integrator('vode', method='bdf', with_jacobian=True)
    r.set_initial_value(init, 0).set_f_params(0.0).set_jac_params(0.0)#.set_solout(solout)

    solout(0, np.ndarray.item(init[states[G]['outputs'][output]].real))

    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        out = states[G]['C']['value']*r.y.reshape(r.y.size,1) + states[G]['D']['value']*inp
        solout(r.t, np.ndarray.item(out[states[G]['outputs'][output]].real))
        print([r.t, np.ndarray.item(out[states[G]['outputs'][output]].real)], file=sys.stderr)

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Reponse", tools=TOOLS,
                   x_axis_label='time', y_axis_label=output, width=800, height=300)

    # add a line renderer with legend and line thickness
    p.line(ts, ys, legend_label=output, line_width=2)

    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states(states)

    return 'Response plot created. Type gdisplay ' + plot + ' to display'

# wrapper method for simulate when finding response initial conditions
@api.dispatcher.add_method
def response(output, G, init, inp, time, plot):
    states = helper.load_states()
    inp = np.matrix(json.loads(inp))
    init = np.matrix(json.loads(init))
    return simulate(output, G, init, inp, time, plot)

# wrapper method for simluate when finding response to a step forcing term
@api.dispatcher.add_method
def step(output, G, control, st, time, plot):
    states = helper.load_states()
    inp = np.zeros((states[G]['B']['value'].shape[1], 1))
    inp[states[G]['controls'][control],0] = st
    init = np.zeros((states[G]['A']['value'].shape[0], 1))
    return simulate(output, G, init, inp, time, plot)
