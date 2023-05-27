from flask import Blueprint
from jsonrpc.backend.flask import api
import numpy as np
import bokeh.plotting as plt
import csv

import helper

data_api = Blueprint('data_api', __name__)
# @simple_page.route('/<page>')
# def show(page):

# states = helper.load_states()

@api.dispatcher.add_method
def data(file, name):
    states[name] = {'value': file, 'meta': {'what': 'data', 'value': file}}
    helper.save_states()
    return 'Data in ' + file + ' will be referenced by: ' + name

@api.dispatcher.add_method
def fit(name, prestatestor):

    x = []
    y = []
    with open(states[name]['value'], mode='r') as infile:
        mappingfile = csv.reader(infile)
        for row in mappingfile:
            x.append(row[0].strip())
            y.append(row[1].strip())

    theta = np.ones(2)
    states[prestatestor] = {'value': theta, 'meta': {'what': 'lin_reg_pred', 'value': theta}}
    helper.save_states()

@api.dispatcher.add_method
def gaussian_process(data_name, plot):

    x = []
    y = []
    with open(states[data_name]['value'], mode='r') as infile:
        mappingfile = csv.reader(infile)
        for row in mappingfile:
            x.append(float(row[0].strip()))
            y.append(float(row[1].strip()))
    xi, xj = np.meshgrid(x, x)

    K = np.exp(-0.5*(xi-xj)*(xi-xj)/0.0001)# + 0.001*np.identity(len(x))

    size = 1000
    #i, j = np.mgrid[:size, :size]/size
    xs = np.linspace(0, 20, size)
    xsi, xsj = np.meshgrid(xs, xs)
    Ks = np.exp(-0.5*(xsi-xsj)*(xsi-xsj)/0.0001)# + 0.001*np.identity(size)

    xs2i, xs2j = np.meshgrid(x, xs)
    Ks2 = np.exp(-0.5*(xs2i-xs2j)*(xs2i-xs2j)/0.0001)

    Kinv = np.linalg.inv(K) #+ 0.001*np.identity(K.shape[0])

    mean = np.asmatrix(Ks2)*np.asmatrix(Kinv)*np.asmatrix(y).transpose()
    cov = np.asmatrix(Ks) - np.asmatrix(Ks2)*Kinv*np.asmatrix(Ks2).transpose()# + 0.01*np.identity(size)

    L = np.linalg.cholesky(cov)
    mu, sigma = 0, 1 # mean and standard deviation
    u = np.random.normal(mu, sigma, (size,1))
    m = np.zeros(size).reshape(size,1)
    ys = mean + np.matrix(L)*np.asmatrix(u)
    ys = np.asarray(ys).reshape(-1)

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Data Visualization", tools=TOOLS,
                   x_axis_label='x', y_axis_label='y')

    # add a line renderer with legend and line thickness
    p.circle(xs,ys, size=1, color="black", alpha=0.9)
    p.square(x,y, size=5, color="blue", alpha=0.2)

    p.quad(top=ys+np.diag(cov), bottom=ys-np.diag(cov), left=xs-0.01,
           right=xs+0.01, color="#B3DE69", alpha = 0.1)

    from bokeh.resources import CDN
    from bokeh.embed import components
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states()

    return 'Data plot created. Type gdisplay ' + plot + ' to display'


@api.dispatcher.add_method
def uniform_gaussian_process(plot):
    size = 500
    i, j = np.mgrid[:size, :size]/size
    x = np.linspace(0, 1, size)
    K = np.exp(-0.5*(i-j)*(i-j)/0.0005) + 0.001*np.identity(size)

    L = np.linalg.cholesky(K)
    mu, sigma = 0, 1 # mean and standard deviation
    u = np.random.normal(mu, sigma, (size,1))
    m = np.zeros(size).reshape(size,1)
    y = m + np.matrix(L)*np.asmatrix(u)
    y = np.asarray(y).reshape(-1)

    print(x, file=sys.stderr)
    print(y, file=sys.stderr)

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Data Visualization", tools=TOOLS,
                   x_axis_label='x', y_axis_label='y')

    # add a line renderer with legend and line thickness
    #p.line(x, y, legend=output, line_width=2)
    p.circle(x,y, size=1, color="navy", alpha=0.5)

    from bokeh.resources import CDN
    from bokeh.embed import components
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states()

    return 'Data plot created. Type gdisplay ' + plot + ' to display'


@api.dispatcher.add_method
def prestatest(data, prestatestor):
    aug_data = np.array([1, data])
    print(aug_data, file=sys.stderr)
    print(states[prestatestor]['value'], file=sys.stderr)
    return np.dot(aug_data, states[prestatestor]['value'])

@api.dispatcher.add_method
def plotdata(name, plot):
    x = []
    y = []
    with open(states[name]['value'], mode='r') as infile:
        mappingfile = csv.reader(infile)
        for row in mappingfile:
            x.append(row[0].strip())
            y.append(row[1].strip())

    # create a new plot with a title and axis labels
    TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"
    p = plt.figure(title="Data Visualization", tools=TOOLS,
                   x_axis_label='x', y_axis_label='y')

    # add a line renderer with legend and line thickness
    #p.line(x, y, legend=output, line_width=2)
    p.circle(x,y, size=10, color="navy", alpha=0.5)

    from bokeh.resources import CDN
    from bokeh.embed import components
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states()

    return 'Data plot created. Type gdisplay ' + plot + ' to display'
