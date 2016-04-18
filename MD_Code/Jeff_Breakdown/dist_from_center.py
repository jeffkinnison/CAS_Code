import numpy as np
import global_variables as gv

def calculate_distance_from_center(center, values):
    """
    Given a ball center, calculate the distance of the supplied trajectory CVs

    Parameters:
        center - an iterable containing the coordinates of the ball center
        values - an iterable containing the Walke coordinates

    Returns:
        The distance between the provided CVs and the ball center as l2 norm
    """
    distance = 0.0
    for i in range(len(center)):
        # standard l2 norm
        if gv.angle_cvs[i] == 0:
            distance += (values[i] - center[i]) ** 2
        else:
            # Adjust dihedrals to fall into [-180,180] deg 
            if values[i] - center[i] > 180.0:
                distance += (values[i] - center[i] - 360.0) ** 2
            elif values[i] - center[i] < -180.0:
                distance += (values[i] - center[i] + 360.0) ** 2
            else:
                distance += (values[i] - center[i]) ** 2
    if abs(distance) < 1.0e-10:
        # Cut off at some min distance
        distance = 0.0
    return np.sqrt(distance)
