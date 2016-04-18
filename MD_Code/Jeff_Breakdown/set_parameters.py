import global_variables as gv

def set_parameters():
    gv.main_directory = p.main_directory
    gv.initial_configuration_directory = p.initial_configuration_directory
    gv.simulation_flag = p.simulation_flag
    gv.balls_flag = p.balls_flag
    gv.rate_flag = p.rate_flag
    gv.num_states = p.num_states
    gv.enhanced_sampling_flag = p.enhanced_sampling_flag
    gv.num_balls_limit = p.num_balls_limit
    gv.radius = p.radius
    gv.num_walkers = p.num_walkers
    gv.num_cvs = p.num_cvs
    gv.lower_bound = p.lower_bound
    gv.upper_bound = p.upper_bound
    gv.angle_cvs = p.angle_cvs
    gv.initial_step_num = p.initial_step_num
    gv.max_num_steps = p.max_num_steps
    gv.num_occupied_balls = p.num_occupied_balls
    gv.first_walker = p.first_walker
    gv.last_walker = p.last_walker
    if gv.enhanced_sampling_flag == 1:
        gv.less_or_greater_flag = p.less_or_greater_flag
        gv.static_threshold_flag = p.static_threshold_flag
        gv.threshold_values = p.threshold_values
        gv.properties_to_keep_track = p.properties_to_keep_track
    elif gv.enhanced_sampling_flag == 2:
        gv.num_occupied_big_clusters = p.num_occupied_big_clusters
        gv.num_occupied_small_clusters = p.num_occupied_small_clusters
        gv.num_balls_for_sc = p.num_balls_for_sc
        gv.num_clusters = p.num_clusters
        gv.num_walkers_for_sc = p.num_walkers_for_sc

    ball_volume = (np.pi ** (gv.num_cvs / 2) * gv.radius ** gv.num_cvs) / special.gamma((gv.num_cvs / 2) + 1)
    if ball_volume != 0.0:
        max_num_balls = int(np.floor((gv.upper_bound - gv.lower_bound) ** gv.num_cvs / ball_volume))*2
    if max_num_balls < gv.num_balls_limit:
        gv.num_balls_limit = max_num_balls
    print 'max # of balls (n_b) = ' + str(gv.num_balls_limit)
    gv.current_num_balls = 0
    if gv.enhanced_sampling_flag == 2 and gv.simulation_flag != 0 and \
                            gv.num_occupied_small_clusters+gv.num_occupied_big_clusters != 0:
        gv.total_num_walkers = gv.num_occupied_small_clusters*gv.num_walkers \
                               + gv.num_occupied_big_clusters*gv.num_walkers_for_sc
    else:
        gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers
    gv.sc_performed = 0