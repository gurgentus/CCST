import os, math, json
import numpy as np
from flask import Blueprint
from jsonrpc.backend.flask import api
import helper
from flask import current_app as app
from PIL import Image
import imageio.v2 as imageio

# orbital mechanics toolbox
import mpc.omt as omt
import nums.OrbitTransfer as OrbitTransfer
import nums.Car as Car

pathplanning_api = Blueprint('pathplanning_api', __name__)

@api.dispatcher.add_method
def generate_grid(name, height, width):
    image_shape = (height, width)
    states = helper.load_states()
    f = os.path.join(app.config['UPLOAD_FOLDER'], states[name]['value'])
    # Load grayscale image and resize using PIL
    img = imageio.imread(f, mode='L')  # 'L' mode for grayscale
    img_pil = Image.fromarray(img)
    img_resized = img_pil.resize((width, height), Image.Resampling.LANCZOS)
    image = np.array(img_resized)
    str_vec = helper.to_str_repr(np.array(image))
    states[name+'_grid'] = {'value': np.array(image), 'meta': {'what': 'matrix', 'value': str_vec}}
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
    states[new_name] = {'value': np.array(array_np), 'meta': {'what': 'matrix', 'value': str_vec}}
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


@api.dispatcher.add_method
def generate_optimal_stop_traj(vel, dist, plot):
    orb_ins = Car.Car()
    orb_ins.SetParams(5,500)
    res = orb_ins.GenerateOptimalStoppingTrajectory([-dist, 0], [vel, 0])
    helper.create_plot(res[0::2], res[1::2], plot, 'Trajectory', 't [s]', 'acceleration control [m/s^2]', 'u(t)')
    return 'Trajectory plot created. Type gdisplay ' + plot + ' to display the trajectory'


@api.dispatcher.add_method
def generate_2d_traj(start_pos, start_vel, start_acc, end_pos, end_vel, end_acc, time_ahead, plot, points=500, fuzzy=False, num=1):
    car_object = Car.Car()
    car_object.SetParams(5,int(points))
    res = car_object.Generate2DTrajectory(json.loads(start_pos), json.loads(start_vel), json.loads(start_acc), json.loads(end_pos), json.loads(end_vel), json.loads(end_acc), float(time_ahead), bool(fuzzy), int(num))
    num_feas_traj = int(len(res)/2);
    if (num_feas_traj < 1):
        return 'No feasible trajectories found.'
    helper.create_multi_plot(res[0::2], res[1::2], plot, 'Trajectory', 's [m]', 'd [m]', '')
    helper.create_plot(res[0], res[1], plot+'_best', 'Best Trajectory', 's [m]', 'd [m]', '')
    return str(num_feas_traj) + ' feasible trajectories generated. Type gdisplay ' + plot + ' to display the trajectories or gdisplay ' + plot + '_best to display the lowest cost trajectory.'


@api.dispatcher.add_method
def orb_transfer_lobatto3(mu, m0, Isp, T, r0, days, timestep_hrs, N, plot):
    orb_ins = OrbitTransfer.OrbitTransfer()
    status = orb_ins.run(float(mu), float(m0), float(Isp), float(T), float(r0), float(days), float(timestep_hrs), int(N))
    if (status == 1):
        return 'Max number of iterations exhausted.'

    if (status == 2):
        return 'Collocation scheme failed.  Try decreasing the timestep.'

    ts = orb_ins.getX()
    ys = orb_ins.getY()
    helper.create_plot(ts, ys, plot, 'Orbit Transfer', 'x [km]', 'y [km]', 'r')
    ts = orb_ins.getT()
    ys = orb_ins.getAngle()
    helper.create_plot(ts, ys, plot, 'Steering Control', 't [sec]', 'angle [deg]', 'steering', height=300)

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

    ts = orb_ins.getX()
    ys = orb_ins.getY()
    helper.create_plot(ts, ys, plot, 'Orbit Transfer', 'x [km]', 'y [km]', 'r')
    ts = orb_ins.getT()
    ys = orb_ins.getAngle()
    helper.create_plot(ts, ys, plot + '_control', 'Steering Control', 't [sec]', 'angle [deg]', 'steering', height=300)


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