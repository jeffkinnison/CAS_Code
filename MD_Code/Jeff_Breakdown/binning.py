def binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers, key_to_ball):
    """Place the walkers into bins.

    Arguments
    """
    initial_weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    initial_weights_array = np.array(initial_weights)
    flux = np.zeros((gv.num_states, gv.num_states))
    flux_num_walkers = np.zeros((gv.num_states, gv.num_states))
    walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights

    start = 0  # indicates whether we are dealing with the very first walker or not

    for i in walker_indices:
        # first, go to walker directory i
        walker_directory = gv.main_directory + '/CAS/walker' + str(i)
        os.chdir(walker_directory)

        # then, obtain new coordinates' values
        if os.path.exists(walker_directory + '/coordinates.out'):
            coordinates = np.loadtxt('coordinates.out')
            if gv.num_cvs > 1:
                new_coordinates = coordinates.tolist()
            else:
                new_coordinates = [float(coordinates)]
            rm_command = 'rm -rf *.out'
            os.system(rm_command)

            # also, write the new coordinates' values  on the trajectory file
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()
        else:
            f = open('trajectory.txt', 'r')
            new_coordinates = f.readlines()[-1].strip().split()
            new_coordinates = [float(coordinate) for coordinate in new_coordinates]
            f.close()

        previous_coordinates = walker_list[i].current_coordinates
        previous_ball_center = walker_list[i].current_ball_center
        previous_distance_from_center = walker_list[i].current_distance_from_center
        initial_step_num = walker_list[i].initial_step_num
        weight = walker_list[i].weight

        if gv.rate_flag == 1:
            state = check_state_function.check_state_function(new_coordinates)
            if walker_list[i].state != -1 and state == -1:
                state = walker_list[i].state
            if walker_list[i].state != -1 and state != -1:
                flux[walker_list[i].state, state] += walker_list[i].weight
                flux_num_walkers[walker_list[i].state, state] += 1
        else:
            state = -1

        inside = 0  # indicates whether we are dealing with the very first walker or not
        # if we're dealing with the very first walker, create the very first ball for the walker
        if (gv.balls_flag == 0 and start == 0) or (gv.balls_flag == 1 and start == 0 and step_num == 0):
            start += 1
            inside += 1
            current_ball_center = [coordinate for coordinate in new_coordinates]
            center_r_key_num = copy.deepcopy(current_ball_center)
            center_r_key_num.append(gv.radius)
            center_r_key_num.append(gv.current_num_balls)
            center_r_key_num.append(1)
            balls[gv.current_num_balls] = np.asarray(center_r_key_num)
            ball_to_walkers[tuple(current_ball_center)] = [i]
            key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
            temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                previous_ball_center, current_ball_center, gv.current_num_balls,
                                                previous_distance_from_center, 0.0, initial_step_num, weight, state)
            gv.current_num_balls += 1

        distance = 0.0
        ball_key = 0
        # otherwise, loop through all of the balls and find the ball that has a center nearest the walker
        if inside == 0:
            for j in range(balls.shape[0]):
                current_ball_center = balls[j][0:gv.num_cvs].tolist()
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                if distance_from_center <= gv.radius or abs(distance_from_center - gv.radius) < 1.0e-10:
                    inside += 1
                if j == 0:
                    distance = distance_from_center
                    ball_key = j
                else:
                    if distance_from_center < distance:
                        distance = distance_from_center
                        ball_key = j

            # walker is inside some ball or is not but needs to be binned to the nearest ball
            if inside != 0 or (inside == 0 and gv.current_num_balls == gv.num_balls_limit):
                balls[ball_key][gv.num_cvs+2] += 1
                current_ball_center = balls[ball_key][0:gv.num_cvs].tolist()
                ball_to_walkers[tuple(current_ball_center)].append(i)
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, ball_key,
                                                    previous_distance_from_center, distance_from_center,
                                                    initial_step_num, weight, state)

            # or walker does not belong in any ball -> create a new ball
            elif inside == 0 and gv.current_num_balls < gv.num_balls_limit:
                current_ball_center = [coordinate for coordinate in new_coordinates]
                center_r_key_num = copy.deepcopy(current_ball_center)
                center_r_key_num.append(gv.radius)
                center_r_key_num.append(gv.current_num_balls)
                center_r_key_num.append(1)
                balls = np.append(balls, [np.asarray(center_r_key_num)], axis=0)
                ball_to_walkers[tuple(current_ball_center)] = [i]
                key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, gv.current_num_balls,
                                                    previous_distance_from_center, 0.0, initial_step_num, weight, state)
                gv.current_num_balls += 1

        # finally, write the new ball on the trajectory file
        current_ball_center = temp_walker_list[i].current_ball_center
        ball_key = temp_walker_list[i].ball_key
        center_r_key_state = copy.deepcopy(current_ball_center)
        center_r_key_state.append(gv.radius)
        center_r_key_state.append(ball_key)
        center_r_key_state.append(state)
        f = open('ball_trajectory.txt', 'a')
        f.write(' '.join(map(lambda coordinate: str(coordinate), center_r_key_state)))
        f.write('\n')
        f.close()

    os.chdir(gv.main_directory + '/CAS')
    np.savetxt('balls_' + str(step_num + 1) + '.txt', balls, fmt=' %+1.5f')
    if gv.rate_flag == 1:
        np.savetxt('flux_' + str(step_num + 1) + '.txt', flux, fmt=' %1.5e')
        np.savetxt('flux_num_walkers_' + str(step_num + 1) + '.txt', flux_num_walkers, fmt=' %d')
    return balls