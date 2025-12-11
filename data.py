"""
Data operations module.
Pure Python functions for data fitting, gaussian processes, and plotting.
"""
import numpy as np
import bokeh.plotting as plt
import csv
import sys
import helper


# Load data from file
def data(file, name):
    """Associate a file with a variable name"""
    states = helper.load_states()
    states[name] = {'value': file, 'meta': {'what': 'data', 'value': file}}
    helper.save_states(states)
    return 'Data in ' + file + ' will be referenced by: ' + name


# Fit linear regression
def fit(name, predictor):
    """Fit a linear regression model to data"""
    states = helper.load_states()

    x = []
    y = []
    with open(states[name]['value'], mode='r') as infile:
        mappingfile = csv.reader(infile)
        for row in mappingfile:
            x.append(row[0].strip())
            y.append(row[1].strip())

    theta = np.ones(2)
    states[predictor] = {'value': theta, 'meta': {'what': 'lin_reg_pred', 'value': theta}}
    helper.save_states(states)


# Gaussian process
def gaussian_process(data_name, plot):
    """Generate gaussian process plot from data"""
    states = helper.load_states()

    x = []
    y = []
    with open(states[data_name]['value'], mode='r') as infile:
        mappingfile = csv.reader(infile)
        for row in mappingfile:
            x.append(float(row[0].strip()))
            y.append(float(row[1].strip()))
    xi, xj = np.meshgrid(x, x)

    K = np.exp(-0.5*(xi-xj)*(xi-xj)/0.0001)

    size = 1000
    xs = np.linspace(0, 20, size)
    xsi, xsj = np.meshgrid(xs, xs)
    Ks = np.exp(-0.5*(xsi-xsj)*(xsi-xsj)/0.0001)

    xs2i, xs2j = np.meshgrid(x, xs)
    Ks2 = np.exp(-0.5*(xs2i-xs2j)*(xs2i-xs2j)/0.0001)

    Kinv = np.linalg.inv(K)

    mean = np.asmatrix(Ks2)*np.asmatrix(Kinv)*np.asmatrix(y).transpose()
    cov = np.asmatrix(Ks) - np.asmatrix(Ks2)*Kinv*np.asmatrix(Ks2).transpose()

    L = np.linalg.cholesky(cov)
    mu, sigma = 0, 1
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

    helper.save_states(states)

    return 'Data plot created. Type gdisplay ' + plot + ' to display'


# Uniform gaussian process
def uniform_gaussian_process(plot):
    """Generate uniform gaussian process plot"""
    states = helper.load_states()

    size = 500
    i, j = np.mgrid[:size, :size]/size
    x = np.linspace(0, 1, size)
    K = np.exp(-0.5*(i-j)*(i-j)/0.0005) + 0.001*np.identity(size)

    L = np.linalg.cholesky(K)
    mu, sigma = 0, 1
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
    p.circle(x,y, size=1, color="navy", alpha=0.5)

    from bokeh.resources import CDN
    from bokeh.embed import components
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states(states)

    return 'Data plot created. Type gdisplay ' + plot + ' to display'


# Prediction test
def predict(data, predictor):
    """Make a prediction using a trained predictor"""
    states = helper.load_states()
    aug_data = np.array([1, data])
    print(aug_data, file=sys.stderr)
    print(states[predictor]['value'], file=sys.stderr)
    return np.dot(aug_data, states[predictor]['value'])


# Plot data from file
def plotdata(name, plot):
    """Plot data from a CSV file"""
    states = helper.load_states()

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
    p.circle(x,y, size=10, color="navy", alpha=0.5)

    from bokeh.resources import CDN
    from bokeh.embed import components
    script, div = components(p)
    states[plot] = {'meta': {}}
    states[plot]['meta'] = {'what': 'plot', 'script': script, 'div': div}

    helper.save_states(states)

    return 'Data plot created. Type gdisplay ' + plot + ' to display'
