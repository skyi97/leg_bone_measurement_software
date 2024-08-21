import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy import optimize
import circle_fit
from scipy.spatial import Voronoi
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, pi


# owned
from msk_gen_fun import SimpleCheck
from msk_vtk_fun import *  # read, extract_points_from_shape, add_append_poly_filter, write, add_mirroring
from msk_vtk import SimpleCoordSystem, SimpleVtkHead, SimpleVtkBody, SimpleVtkTail, \
    ComplexVtkExtraction, SimpleVtkExtraction, \
    SimpleVtkDefaultShape
from msk_align import SimpleAlignSVDFemur, ComplexAlignSVDLowerlimb
from msk_math import *
from _02_femur_functions import*

def least_squares_circle(coords): #taken from circle_fit library
    """not used"""
    """
    Circle fit using least-squares solver.
    Inputs:

        - coords, list or numpy array with len>2 of the form:
        [
    [x_coord, y_coord],
    ...,
    [x_coord, y_coord]
    ]
        or numpy array of shape (n, 2)

    Outputs:

        - xc : x-coordinate of solution center (float)
        - yc : y-coordinate of solution center (float)
        - R : Radius of solution (float)
        - residu : MSE of solution against training data (float)
        *residu has been modified by Sky to give standard deviation instead
    """

    x, y = None, None
    if isinstance(coords, np.ndarray):
        x = coords[:, 0]
        y = coords[:, 1]
    elif isinstance(coords, list):
        x = np.array([point[0] for point in coords])
        y = np.array([point[1] for point in coords])
    else:
        raise Exception("Parameter 'coords' is an unsupported type: " + str(type(coords)))

    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, _ = optimize.leastsq(f, center_estimate, args=(x,y))
    xc, yc = center
    Ri = calc_R(x, y, *center)
    R = Ri.mean()

    residu = np.mean(Ri - R)

    return xc, yc, R, residu




def maximum_inscribed_circle(sectioned_points, point_whole_shape):
    """not used"""
    """uses voronoi method to calculate the midpoint and radius of the smallest circle that
    can be fit within points"""
    def closest_node(node, nodes):
        nodes = np.asarray(nodes)
        deltas = nodes - node
        dist = np.linalg.norm(deltas, axis=1)
        min_idx = np.argmin(dist)
        return nodes[min_idx], dist[min_idx], deltas[min_idx][1] / deltas[min_idx][0]  # point, distance, slope

    points = sectioned_points
    plane_20mm_numpy_2D = point_whole_shape

    vor = Voronoi(points)

    max_d = 0
    max_v = None
    vertices = vor.vertices
    vertices = vertices[vertices[:, 0] < max(plane_20mm_numpy_2D[:, 0])]
    vertices = vertices[min(plane_20mm_numpy_2D[:, 0]) < vertices[:, 0]]

    vertices = vertices[vertices[:, 1] < max(plane_20mm_numpy_2D[:, 1])]
    vertices = vertices[min(plane_20mm_numpy_2D[:, 1]) < vertices[:, 1]]

    for v in vertices:

        _, d, _ = closest_node(v, points)
        if d > max_d:
            max_d = d
            max_v = v

    xc = max_v[0]
    yc = max_v[1]
    r = max_d
    sd = np.std(np.sqrt((sectioned_points[:, 0] - xc) ** 2 + (sectioned_points[:, 1] - yc) ** 2))

    return xc, yc, r, sd


def furthest_points_in_slice(points):  # measure distance between every points in the plane and find the max distance

    arr1 = points
    arr2 = points

    max_diameter = [[0, 0, 0], [0, 0, 0], 0]

    for x in arr1:
        for y in arr2:
            dist = np.sqrt(np.sum((x - y) ** 2, axis=0))
            if dist > max_diameter[2]:
                max_diameter[0] = x
                max_diameter[1] = y
                max_diameter[2] = dist
    print(max_diameter)
    center = np.add(max_diameter[0], max_diameter[1]) / 2

    # radius = np.sqrt(s_in[:,0]**2 + s_in[:,1]**2 - 2*s_in[:,0]*x_center - 2*s_in[:,1]* y_center + x_center**2 + y_center **2)

    return max_diameter


def mediolateral_dimensions(s_in):
    proximal = SimpleVtkExtraction(s_in).get_cube_by_range([0, 100])
    distal = SimpleVtkExtraction(s_in).get_cube_by_range([100, 0])

    distal_numpy = extract_points_from_shape(distal)
    proximal_numpy = extract_points_from_shape(proximal)

    dist_lateral = distal_numpy[distal_numpy[:, 0] == max(distal_numpy[:, 0])].flatten()
    dist_medial = distal_numpy[distal_numpy[:, 0] == min(distal_numpy[:, 0])].flatten()

    dist_ML_midpoint = (dist_lateral + dist_medial) / 2

    prox_lateral = proximal_numpy[proximal_numpy[:, 0] == max(proximal_numpy[:, 0])].flatten()
    prox_medial = proximal_numpy[proximal_numpy[:, 0] == min(proximal_numpy[:, 0])].flatten()

    prox_ML_midpoint = (prox_lateral + prox_medial) / 2

    #### visualisation section
    # dist_ML = SimpleVtkDefaultShape().get_cylinder_by_points(1, dist_lateral, dist_medial)
    # prox_ML = SimpleVtkDefaultShape().get_cylinder_by_points(1, prox_lateral, prox_medial)
    # tibia_aa = SimpleVtkDefaultShape().get_cylinder_by_points(1, dist_ML_midpoint, prox_ML_midpoint)
    # visualise([s_in, dist_ML, prox_ML, tibia_aa])

    #########

    return prox_ML_midpoint, dist_ML_midpoint, prox_lateral, prox_medial, dist_lateral, dist_medial


def rotate_bone(s_in, v1, v2):
    """not used"""
    """ let v1 be the target vector to be rotated to e.g anatomical axis
    v2 should be the desired coordinate to be rotated to, e.g z axis"""


    uv1 = v1
    uv2 = v2


    x_disp = uv1[0] - uv2[0]
    y_disp = uv1[1] - uv2[1]
    z_disp = uv1[2] - uv2[2]
    lis = [x_disp,y_disp,z_disp]
    print(lis)
    new_lis = []
    for i in lis:


        if i == 0:
            new_lis.append(0)
        if i< 0:
            new_lis.append(-1)
        else:
            new_lis.append(1)



    threeD, coronal, saggital, transverse = cal_angle_by_vectors(v1, v2)
    print(new_lis)
    coronal = coronal*new_lis[0]
    saggital = saggital* new_lis[1]
    transverse = transverse * new_lis[2]

    rotate_angle = [transverse, -coronal, saggital]

    # rotate_angle = cal_angel_by_points(v1,[0,0,0],v2,[0,0,0],angle_tpye = 'eulZYX')
    print('rotate angle: ', rotate_angle)
    # rotate_angle[0] -= 90
    # rotate_angle = -rotate_angle

    s_in_rotated = add_transform_by_rotation(s_in, rotation=rotate_angle)

    return s_in_rotated



def get_height(points):  # returns z to z distance of most proximal and most distal point
    # s_in = SimpleAlignSVDFemur(f_in, isRight=True).get_aligned_femur()
    # points = extract_points_from_shape(s_in)
    point_proxima = points[points[:, 2] == max(points[:, 2])].flatten()  # points[:,2] accesses the column index 2
    point_distal = points[points[:, 2] == min(points[:, 2])].flatten()
    z_length = point_proxima[2] - point_distal[2]

    return [z_length, point_proxima, point_distal]

