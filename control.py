"""
Control systems operations module.
Pure Python functions for matrix operations, state-space systems, and control theory.
"""
import numpy as np
import bokeh.plotting as plt
import csv, json, sys
from scipy import linalg as la
from scipy.integrate import ode
from bokeh.resources import CDN
from bokeh.embed import components
import helper


# Define a matrix
def matrix(name, value):
    """Create and store a matrix by name"""
    states = helper.load_states()
    val = np.array(json.loads(value) if isinstance(value, str) else value)
    str_vec = helper.to_str_repr(val)

    states[name] = {'value': val, 'meta': {'what': 'matrix', 'value': str_vec}}
    helper.save_states(states)

    return states[name]['value'].tolist()


# Calculate eigenvectors/eigenvalues of a matrix
def modes(mat, V, E):
    """
    Calculate eigenvalues and eigenvectors of a matrix.

    Args:
        mat: Name of matrix in states or matrix value directly
        V: Name to store eigenvectors
        E: Name to store eigenvalues

    Returns:
        Eigenvalues as string representation
    """
    states = helper.load_states()

    # Support both matrix name and direct matrix value
    if isinstance(mat, str):
        matrix_val = states[mat]['value']
    else:
        matrix_val = mat

    e_vals, e_vecs = np.linalg.eig(matrix_val)
    str_vec = helper.to_str_repr(e_vecs)
    str_vals = helper.to_str_repr(e_vals)
    states[V] = {'value': np.array(e_vecs), 'meta': {'what': 'matrix', 'value': str_vec}}
    states[E] = {'value': e_vals, 'meta': {'what': 'vector', 'value': str_vals}}

    helper.save_states(states)

    return states[E]['meta']['value']


# Define state-space system using matrices A,B,C,D and store it under variable name
def ss(A, B, C, D, name):
    """Create a state-space system"""
    states = helper.load_states()
    states[name] = {'value': [A, B, C, D], 'A': states[A], 'B': states[B], 'C': states[C],
     'D': states[D], 'controls': {}, 'outputs': {}, 'meta': {'what': 'ss',
     'A': states[A]['meta'], 'B': states[B]['meta'], 'C': states[C]['meta'], 'D': states[D]['meta']}}
    helper.save_states(states)
    result = 'x\' = ' + A + 'x + ' + B + 'u , y = ' + C + 'x + ' + D + 'u'
    return result


# Define a controller using matrices A, B1, B2, C, D1, D2 and store under name
def controller(A, B1, B2, C, D1, D2, name):
    """Create a controller system"""
    states = helper.load_states()
    states[name] = {'value': [A, B1, B2, C, D1, D2], 'A': states[A], 'B1': states[B1], 'B2': states[B2],
     'C': states[C], 'D1': states[D1], 'D2': states[D2], 'controls': {}, 'outputs': {},
     'meta': {'what': 'controller',
     'A': states[A]['meta'], 'B1': states[B1]['meta'], 'B2': states[B2]['meta'], 'C': states[C]['meta'],
     'D1': states[D1]['meta'], 'D2': states[D2]['meta']}}
    helper.save_states(states)
    result = 'x_c\' = ' + A + 'x_c + ' + B1 + 'y + ' + B2 + 'r , u = ' + C + 'x_c + ' + D1 + 'y + ' + D2 + 'r'
    return result


# Connect state space system G and a controller K into a feedback loop and store
# the resulting system under name
def feedback(G, K, name):
    """Create a feedback control system"""
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


# Associate the num's control variable with name
def control(num, G, name):
    """Set control variable for a system"""
    states = helper.load_states()
    states[G]['controls'][name] = num - 1
    helper.save_states(states)
    return 'Control number ' + str(num) + ' of ' + G + ' is set to ' + name


# Associate the num's output variable with name
def output(num, G, name):
    """Set output variable for a system"""
    states = helper.load_states()
    states[G]['outputs'][name] = num - 1
    helper.save_states(states)
    return 'Output number ' + str(num) + ' of ' + G + ' is set to ' + name


# Simulate output response of the state-space system G to initial conditions init
# and a constant forcing term inp for a time period time and save the
# resulting graph to plot
def simulate(output, G, init, inp, time, plot):
    """Simulate system response"""
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
    r.set_initial_value(init, 0).set_f_params(0.0).set_jac_params(0.0)

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


# Wrapper method for simulate when finding response initial conditions
def response(output, G, init, inp, time, plot):
    """Calculate system response to initial conditions and input"""
    states = helper.load_states()
    inp = np.array(json.loads(inp)) if isinstance(inp, str) else inp
    init = np.array(json.loads(init)) if isinstance(init, str) else init
    return simulate(output, G, init, inp, time, plot)


# Wrapper method for simulate when finding response to a step forcing term
def step(output, G, control_name, st, time, plot):
    """Calculate system step response"""
    states = helper.load_states()
    inp = np.zeros((states[G]['B']['value'].shape[1], 1))
    inp[states[G]['controls'][control_name],0] = st
    init = np.zeros((states[G]['A']['value'].shape[0], 1))
    return simulate(output, G, init, inp, time, plot)
