import global_variables as gv
import walker

def initialize(input_initial_values_file, walker_list, temp_walker_list, balls, ball_to_walkers, vacant_walker_indices):
    """
    Initialize the system to some provided state. 
    
    Parameters:
        input_initial_values_file - file containing the initial values (CVs?) for the walkers
        walker_list - an iterable containing the Walker objects to be used
        temp_walker_list - ?
        balls - the current balls (states) in the system
        ball_to_walkers- ?
        vacant_walker_indices - list of Walkers with no state?
    """
    # Essentially "zero-out" the Walkers'CVs
    for i in range(len(walker_list)):
        walker_list[i] = walker.Walker([-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, i, 0.0, [-1000.0] * gv.num_cvs,
                                       [-1000.0] * gv.num_cvs, 0, 0.0, 0.0, 0, 0.0, -1)
    
    # Handle starting a run from scratch
    if gv.simulation_flag == 0:  # new simulation
        # Normalize all Walker weights to be equal
        initial_weight = 1.0/gv.total_num_walkers

        # Read in the the CVs for eah Walker in the system
        f = open(input_initial_values_file, 'r')
        for n in range(gv.num_occupied_balls):
            initial_values = [None] * gv.num_cvs
            for i in range(gv.num_cvs):
                initial_values[i] = float(f.readline())
            
            # If the macro-state assignment is not known, get it
            if gv.rate_flag == 1:
                initial_state = check_state_function.check_state_function(initial_values)
                print initial_state

            # ???
            for i in range(n * gv.num_walkers, (n + 1) * gv.num_walkers):
                walker_list[i].set(initial_values, initial_weight)
                if gv.rate_flag == 1:
                    walker_list[i].state = initial_state
        f.close()

        # Set up the directory structure and move Walker info to the correct place
        os.system('mkdir CAS')
        os.chdir(gv.main_directory + '/CAS')
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            shutil.copytree(gv.initial_configuration_directory, walker_directory)

    # Handle starting the run in the case that simulations were still running
    elif gv.simulation_flag == 1:  # restarting simulation in the middle of simulation
        # Get the most recent state for each Walker
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            trajectory = np.loadtxt('trajectory.txt')
            if gv.num_cvs == 1:
                previous_coordinates = [trajectory[-2]]
                current_coordinates = [trajectory[-1]]
            else:
                previous_coordinates = trajectory[-2].tolist()
                current_coordinates = trajectory[-1].tolist()
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            ball_trajectory = np.loadtxt('ball_trajectory.txt')
            previous_ball_center = ball_trajectory[-2][0:gv.num_cvs].tolist()
            current_ball_center = ball_trajectory[-1][0:gv.num_cvs].tolist()
            current_ball_radius = ball_trajectory[-1][gv.num_cvs]
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].radius = current_ball_radius
            walker_list[i].previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                          previous_ball_center)
            walker_list[i].current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                          current_ball_center)
            # Zero-out the balls to calculate transition rates
            if gv.rate_flag == 1:
                walker_list[i].state = int(ball_trajectory[-1][-1])

            # Retrieve previous ball definitions if they are directed to persist
            if gv.balls_flag == 1:
                os.chdir(gv.main_directory + '/CAS')
                balls = np.loadtxt('balls_' + str(gv.initial_step_num) + '.txt')
                gv.current_num_balls = balls.shape[0]
    
    # Set up for the case that the run is starting while determining Walker-ball assignments
    elif gv.simulation_flag == 2:  # restarting simulation in the middle of binning
        # Load most recent Walker states
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            trajectory = np.loadtxt('trajectory.txt')
            previous_coordinates = trajectory[-2].tolist()
            current_coordinates = trajectory[-1].tolist()
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            num_lines = sum(1 for line in open('ball_trajectory.txt'))
            # if walker is already binned to a ball, delete the binning and have binning start from scratch
            if num_lines > gv.initial_step_num:
                os.system('sed -i \'$d\' ball_trajectory.txt')
            ball_trajectory = np.loadtxt('ball_trajectory.txt')
            previous_ball_center = ball_trajectory[-2][0:gv.num_cvs].tolist()
            current_ball_center = ball_trajectory[-1][0:gv.num_cvs].tolist()
            current_ball_radius = ball_trajectory[-1][gv.num_cvs]
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].radius = current_ball_radius
            walker_list[i].previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                          previous_ball_center)
            walker_list[i].current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                         current_ball_center)
            # Prep the walkers for transition rate calculations
            if gv.rate_flag == 1:
                walker_list[i].state = int(ball_trajectory[-1][-1])
            
            # Recover the persistent ball definitions
            if gv.balls_flag == 1:
                os.chdir(gv.main_directory + '/CAS')
                balls = np.loadtxt('balls_' + str(gv.initial_step_num) + '.txt')
                gv.current_num_balls = balls.shape[0]

    # Set up for the case that the run was perviously performing resampling
    elif gv.simulation_flag == 3:  # restarting simulation in the middle of resampling
        # Zero out the total weight for later recovery 
        total_weight = 0.0
        previous_ball_to_walkers = {}
        os.chdir(gv.main_directory + '/CAS')
        previous_balls_weights = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')
        previous_balls_walker_count = np.zeros((previous_balls_weights.shape[0], previous_balls_weights.shape[1]))
        for i in range(previous_balls_weights.shape[0]):
            previous_balls_walker_count[i] = previous_balls_weights[i]
            previous_balls_walker_count[i][-1] = gv.num_walkers

        # TODO: make sure that gv.num_occupied_balls is equal to the highest walker number inside the CAS folder

        for i in range(gv.num_occupied_balls + 1):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            # if all of the files exist in the walker folder, it is a complete walker
            if os.path.isfile(walker_directory + '/weight_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/ball_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/trajectory.txt') and \
                    os.path.isfile(walker_directory + '/traj.xtc') and \
                    os.path.isfile(walker_directory + '/minim.gro') and \
                    os.path.isfile(walker_directory + '/minim.tpr'):
                os.chdir(walker_directory)
                f = open('weight_trajectory.txt', 'r')
                weight = float(f.readlines()[-1].strip())
                walker_list[i].weight = weight
                total_weight += weight
                f.close()

                ball_trajectory = np.loadtxt('ball_trajectory.txt')
                previous_ball = ball_trajectory[-2].tolist()
                previous_ball_key = previous_ball[gv.num_cvs+1]
                previous_ball_center = previous_ball[0:gv.num_cvs]
                previous_balls_weights[previous_ball_key][-1] -= weight
                if previous_balls_weights[previous_ball_key][-1] < 0.0:
                    print 'ERROR: weight is ' + str(previous_balls_weights[previous_ball_key][-1]) + ' for walker ' + \
                          str(i) + ' with ball_key ' + str(previous_ball_key)
                previous_balls_walker_count[previous_ball_key][-1] -= 1
                if previous_balls_walker_count[previous_ball_key][-1] < 0:
                    print 'ERROR: walker count is ' + str(previous_balls_walker_count[previous_ball_key][-1]) + \
                          ' for walker ' + str(i) + ' with ball key ' + str(previous_ball_key)
                if tuple(previous_ball_center) in previous_ball_to_walkers:
                    previous_ball_to_walkers[tuple(previous_ball_center)].append(i)
                else:
                    previous_ball_to_walkers[tuple(previous_ball_center)] = [i]
                current_ball = ball_trajectory[-1].tolist()
                current_ball_radius = current_ball[gv.num_cvs]
                current_ball_key = current_ball[gv.num_cvs+1]
                current_ball_center = current_ball[0:gv.num_cvs]
                if tuple(current_ball_center) in ball_to_walkers:
                    ball_to_walkers[tuple(current_ball_center)].append(i)
                else:
                    ball_to_walkers[tuple(current_ball_center)] = [i]
                    gv.current_num_balls += 1
                walker_list[i].radius = current_ball_radius
                walker_list[i].previous_ball_center = previous_ball_center
                walker_list[i].current_ball_center = current_ball_center

                trajectory = np.loadtxt('trajectory.txt')
                previous_coordinates = trajectory[-2].tolist()
                current_coordinates = trajectory[-1].tolist()
                walker_list[i].previous_coordinates = previous_coordinates
                walker_list[i].current_coordinates = current_coordinates
                if gv.rate_flag == 1:
                    current_state = int(current_ball[-1])
                else:
                    current_state = -1
                walker_list[i].state = current_state
                walker_list[i].ball_key = current_ball_key
                previous_distance_from_center = calculate_distance_from_center(previous_coordinates, previous_ball_center)
                current_distance_from_center = calculate_distance_from_center(current_coordinates, current_ball_center)
                walker_list[i].previous_distance_from_center = previous_distance_from_center
                walker_list[i].current_distance_from_center = current_distance_from_center

                temp_walker_list[i] = walker.Walker(previous_coordinates, current_coordinates, i, current_ball_radius,
                                                    previous_ball_center, current_ball_center, current_ball_key,
                                                    previous_distance_from_center, current_distance_from_center, 0,
                                                    weight, current_state)

            # otherwise, it is an incomplete walker that needs missing files
            else:
                if os.path.isdir(walker_directory):
                    os.chdir(gv.main_directory + '/CAS')
                    os.system('rm -rf walker' + str(i))
                vacant_walker_indices.append(i)

        # create new walkers for the remaining weights
        excess_index = gv.num_occupied_balls + 1
        for i in range(previous_balls_weights.shape[0]):
            if previous_balls_weights[i][-1] > 0.0:
                if previous_balls_walker_count[i][-1] <= 0:
                    print 'ERROR: at least one walker should exist if there is a weight of ' + \
                          str(previous_balls_weights[i][-1]) + ' for walker ' + str(i)
                else:
                    current_ball_center = previous_balls_weights[i][0:gv.num_cvs].tolist()
                    reference_walker = ball_to_walkers[tuple(current_ball_center)][0]
                    reference_walker_directory = gv.main_directory + '/CAS/walker/' + str(reference_walker)
                    if len(vacant_walker_indices) > 0:
                        walker_index = vacant_walker_indices.pop(0)
                    else:
                        walker_index = excess_index
                        excess_index += 1
                    walker_directory = gv.main_directory + '/CAS/walker' + str(walker_index)
                    shutil.copytree(reference_walker_directory, walker_directory)

                    weight = previous_balls_weights[i][-1]
                    previous_balls_weights[i][-1] -= weight

                    os.chdir(walker_directory)
                    f = open('weight_trajectory.txt', 'w')
                    f.write(str(weight) + '\n')
                    walker_list[walker_index].weight = weight
                    total_weight += weight
                    f.close()

                    ball_to_walkers[tuple(current_ball_center)].append(walker_index)
                    walker_list[walker_index].current_ball_center = current_ball_center

                    trajectory = np.loadtxt('trajectory.txt')
                    previous_coordinates = trajectory[-2].tolist()
                    current_coordinates = trajectory[-1].tolist()
                    walker_list[walker_index].previous_coordinates = previous_coordinates
                    walker_list[walker_index].current_coordinates = current_coordinates
                    ball_trajectory = np.loadtxt('ball_trajectory.txt')
                    previous_ball_center = ball_trajectory[-2][0:gv.num_cvs].tolist()
                    walker_list[i].previous_ball_center = previous_ball_center
                    current_state = ball_trajectory[-1][-1]
                    current_ball_key = ball_trajectory[-1][gv.num_cvs+1]
                    current_ball_radius = ball_trajectory[-1][gv.num_cvs]
                    walker_list[walker_index].state = current_state
                    walker_list[walker_index].ball_key = current_ball_key
                    walker_list[walker_index].radius = current_ball_radius
                    previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                       previous_ball_center)
                    current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                      current_ball_center)
                    walker_list[i].previous_distance_from_center = previous_distance_from_center
                    walker_list[i].current_distance_from_center = current_distance_from_center

                    temp_walker_list[walker_index] = walker.Walker(previous_coordinates, current_coordinates,
                                                                   walker_index, current_ball_radius,
                                                                   previous_ball_center, current_ball_center,
                                                                   current_ball_key, previous_distance_from_center,
                                                                   current_distance_from_center, 0, weight,
                                                                   current_state)

        # check if total weight is 1.0
        if total_weight != 1.0:
            print 'ERROR: total weight is ' + str(total_weight)
        if gv.balls_flag == 1:
            os.chdir(gv.main_directory + '/CAS')
            balls = np.loadtxt('balls_' + str(gv.initial_step_num) + '.txt')
            gv.current_num_balls = balls.shape[0]
    return balls
