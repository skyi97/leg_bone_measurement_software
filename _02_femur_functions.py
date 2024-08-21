import numpy as np

# owned

from msk_vtk_fun import *  # read, extract_points_from_shape, add_append_poly_filter, write, add_mirroring
from msk_vtk import SimpleCoordSystem, SimpleVtkHead, SimpleVtkBody, SimpleVtkTail, \
    ComplexVtkExtraction, SimpleVtkExtraction, \
    SimpleVtkDefaultShape
from msk_align import SimpleAlignSVDFemur, ComplexAlignSVDLowerlimb
from msk_math import *  # cal_svd_of_shape, cal_fit_sphere, cal_angle_by_vectors, create_fit_circle



def cal_fna(s_proxima):
    points_femur_neck, _, _ = cal_svd_of_shape(s_proxima, sampling=4)
    v_femur_neck = points_femur_neck[1] - points_femur_neck[0]
    l_femural_neck = SimpleVtkDefaultShape().get_cylinder_by_points(2, points_femur_neck[1], points_femur_neck[0])
    return l_femural_neck, v_femur_neck


def cal_faa(s_shaft):
    points_mid_shaft, _, _ = cal_svd_of_shape(s_shaft, sampling=4)
    print(points_mid_shaft)
    v_f_anatomic = points_mid_shaft[1] - points_mid_shaft[0]
    l_anatomic = SimpleVtkDefaultShape().get_cylinder_by_points(
        1, points_mid_shaft[0]-2*v_f_anatomic, points_mid_shaft[-1]+v_f_anatomic
    )
    return l_anatomic, v_f_anatomic


def cal_fhc(s_fh):
    centre, radius = cal_fit_sphere(extract_points_from_shape(s_fh), plot_check=False)
    s_fh_fit = SimpleVtkDefaultShape().get_sphere(centre, radius)
    l_mechanical = SimpleVtkDefaultShape().get_cylinder_by_points(1, centre, [0, 0, 0])
    return l_mechanical, s_fh_fit, radius,centre


def cal_fst(s_shaft_top):
    points_femur_shaft_top, _, _ = cal_svd_of_shape(s_shaft_top, sampling=4)
    v_femur_shaft_top = points_femur_shaft_top[1] - points_femur_shaft_top[0]
    l_femural_shaft_top = SimpleVtkDefaultShape().get_cylinder_by_points(1, points_femur_shaft_top[1]+0.7*v_femur_shaft_top,
                                                                         points_femur_shaft_top[0])
    return l_femural_shaft_top, v_femur_shaft_top


def cal_fsb(s_shaft_bot):
    points_femur_shaft_bot, _, _ = cal_svd_of_shape(s_shaft_bot, sampling=4)
    v_femur_shaft_bot = points_femur_shaft_bot[1] - points_femur_shaft_bot[0]
    # l_femural_shaft_bot = SimpleVtkDefaultShape().get_line(points_femur_shaft_bot[1], points_femur_shaft_bot[0])
    l_femural_shaft_bot = SimpleVtkDefaultShape().get_cylinder_by_points(1, points_femur_shaft_bot[1],
                                                                         points_femur_shaft_bot[0]-0.7*v_femur_shaft_bot)
    return l_femural_shaft_bot, v_femur_shaft_bot


def cal_diaphyseal_condyle(s_condyle):
    points = extract_points_from_shape(s_condyle)  # vtk to point cloud
    # this code based on a right femur
    points_medial = points[points[:, 0] > 0, :]  # extract the point with x coordinate >0, grouped as medial point
    points_lateral = points[points[:, 0] < 0, :]  # so <0 are lateral points
    points_medial_condyle = points_medial[points_medial[:, 2] == min(points_medial[:, 2]),
                            :].flatten()  # the point with minumum z-coordinate is the points for condyle (most distal point)
    points_lateral_condyle = points_lateral[points_lateral[:, 2] == min(points_lateral[:, 2]), :].flatten()
    # l_diaphyseal_condyle = SimpleVtkDefaultShape().get_line(points_medial_condyle, points_lateral_condyle)
    l_diaphyseal_condyle = SimpleVtkDefaultShape().get_cylinder_by_points(2, points_medial_condyle,
                                                                          points_lateral_condyle)
    v_disphyseal_condyle = points_lateral_condyle - points_medial_condyle
    return l_diaphyseal_condyle, v_disphyseal_condyle


def cal_posterior_condyle(s_condyle):
    points = extract_points_from_shape(s_condyle)
    points_medial = points[points[:, 0] > 0, :]
    points_lateral = points[points[:, 0] < 0, :]
    points_medial_posterior = points_medial[points_medial[:, 1] == max(points_medial[:, 1]), :].flatten()
    points_lateral_posterior = points_lateral[points_lateral[:, 1] == max(points_lateral[:, 1]), :].flatten()

    l_posterior_condyle = SimpleVtkDefaultShape().get_cylinder_by_points(2, points_medial_posterior,
                                                                         points_lateral_posterior)
    v_posterior_condyle = points_medial_posterior - points_lateral_posterior
    return l_posterior_condyle, v_posterior_condyle


def cal_meidaleteral_langth(s_condyle):
    points = extract_points_from_shape(s_condyle)  # vtk to point cloud
    # this code based on a right femur
    points_medial = points[points[:, 0] > 0, :]
    points_lateral = points[points[:, 0] < 0, :]
    points_medial_condyle = points_medial[points_medial[:, 2] == min(points_medial[:, 2]),
                            :].flatten()
    points_lateral_condyle = points_lateral[points_lateral[:, 2] == min(points_lateral[:, 2]), :].flatten()

    squared_dist = np.sum((points_medial_condyle - points_lateral_condyle) ** 2, axis=0)
    dist = np.sqrt(squared_dist)
    return dist


def cal_femur_length(s_bone):
    points = extract_points_from_shape(s_bone)  # vtk to point cloud
    # this code based on a right femur
    point_proxima = points[points[:, 2] == max(points[:, 2])].flatten()
    point_distal = points[points[:, 2] == min(points[:, 2])].flatten()
    print(point_distal)
    squared_dist = np.sum((point_proxima - point_distal) ** 2, axis=0)
    dist = np.sqrt(squared_dist)
    return dist