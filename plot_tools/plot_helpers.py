import numpy as np


def convert_cartesian(position_data: tuple, model_type=None):
    """
    A helper function to convert the coordinates from Spherical or Cylindrical into Cartesian
    The 'q' name scheme is the general 3-D coordinate, they associate with the same convention for what type
    is selected. (E.g. type = sphere, then q1 = radius, q2 = theta and q3 = phi)
    :param position_data:
    :param model_type:
    :return:
    """
    q1 = position_data[0]
    q2 = position_data[1]
    q3 = position_data[2]
    if model_type is None:
        raise OSError("Define co-ordinate system to transform from")
    if model_type == 'is_sphere':
        x = q1 * np.cos(q3) * np.cos(q2)
        y = q1 * np.sin(q3) * np.cos(q2)
        z = q1 * np.sin(q2)
    if model_type == 'is_cylindrical':
        x = q1 * np.cos(q2)
        y = q1 * np.sin(q2)
        z = q3
    return x, y, z


