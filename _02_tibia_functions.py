import numpy as np
import circle_fit

# owned
from msk_vtk_fun import*
from msk_math import*  #
from _03_misc_functions import*




def tibia_shaft_diameters(s_in, points, height):
    std_length = height[0]
    point_proximal = height[1]
    point_distal = height[2]

    plane_80 = SimpleVtkExtraction(s_in).get_plane_view([0, 0, point_distal[2] + std_length * 0.2], [0, 0, 1])
    plane_50 = SimpleVtkExtraction(s_in).get_plane_view([0, 0, point_distal[2] + std_length * 0.5], [0, 0, 1])
    plane_20 = SimpleVtkExtraction(s_in).get_plane_view([0, 0, point_distal[2] + std_length * 0.8], [0, 0, 1])

    plane_80_num = extract_points_from_shape(plane_80)
    plane_50_num = extract_points_from_shape(plane_50)
    plane_20_num = extract_points_from_shape(plane_20)

    len_dict_80 = furthest_points_in_slice(plane_80_num)
    len_dict_50 = furthest_points_in_slice(plane_50_num)
    len_dict_20 = furthest_points_in_slice(plane_20_num)

    return len_dict_20, len_dict_50, len_dict_80


def taa_svd(s_in):
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_in.GetBounds()
    height = zmax - zmin

    s_shaft = SimpleVtkExtraction(s_in).get_cube_by_range([70, height - 120], inside_out=True)
    points_mid_shaft, x_axis, _ = cal_svd_of_shape(s_shaft, sampling=4)

    v_f_anatomic = points_mid_shaft[1] - points_mid_shaft[0]

    top_shaft = SimpleVtkExtraction(s_in).get_plane_view([0, 0, zmax - 100], [0, 0, 1])
    bot_shaft = SimpleVtkExtraction(s_in).get_plane_view([0, 0, zmin + 100], [0, 0, 1])
    n_top = extract_points_from_shape(top_shaft)
    n_bot = extract_points_from_shape(bot_shaft)
    n_top = np.mean(n_top, axis=0).flatten()
    n_bot = np.mean(n_bot, axis=0).flatten()
    v_f_anatomic = n_top - n_bot
    points_mid_shaft = [n_top, n_bot]
    s_shaft = SimpleVtkExtraction(s_in).get_cube_by_range([100, 100], inside_out=True)

    return v_f_anatomic, points_mid_shaft, s_shaft


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

    return prox_ML_midpoint, dist_ML_midpoint, prox_lateral, prox_medial, dist_lateral, dist_medial


def tibia_plateau_axis(s_in, height, screenshot=False):

    plane_20mm = SimpleVtkExtraction(s_in).get_plane_view([0, 0, height[1][2] - 20], [0, 0, 1])

    mcc_20, lcc_20 = condylar_centre_best_fit(plane_20mm)

    vector = lcc_20[0] - mcc_20[0]

    vector[2] = 0
    midpoint = (lcc_20[0] + mcc_20[0]) / 2
    cc_r = [lcc_20[0], lcc_20[1], mcc_20[0], mcc_20[1]]

    v1 = [vector, midpoint, cc_r]
    cc = [lcc_20, mcc_20]

    if screenshot:
        new_axis = SimpleVtkDefaultShape().get_cylinder_by_points(1, mcc_20[0], lcc_20[0])
        plateau = SimpleVtkExtraction(s_in).get_cube_by_range([0, 30])
        lcc_centroid = SimpleVtkDefaultShape().get_cylinder_by_points(lcc_15[1], lcc_20[0],
                                                                      np.array(lcc_20[0]) + [0, 0, 0.01])
        mcc_centroid = SimpleVtkDefaultShape().get_cylinder_by_points(mcc_20[1], mcc_20[0],
                                                                      np.array(mcc_20[0]) + [0, 0, 0.01])

        vtk_input = [plateau, plane_15mm, plane_20mm, new_axis, mcc_centroid, lcc_centroid]

        screenshots(vtk_input, mode='all', autoscale=True)

    return v1


def condylar_centre_best_fit(slice):
    n_slice = extract_points_from_shape(slice)
    n_slice_2D = n_slice[:, [0, 1]]
    z_value = n_slice[0, 2]

    slice_width = max(n_slice_2D[:, 0]) - min(n_slice_2D[:, 0])
    slice_height = max(n_slice_2D[:, 1]) - min(n_slice_2D[:, 1])

    lat_max_width = n_slice_2D[n_slice_2D[:, 0] >
                               (max(n_slice_2D[:, 0]) - 0.4 * slice_width)]
    lat_max_height = max(lat_max_width[:, 1]) - min(lat_max_width[:, 1])
    lat_max_width = max(lat_max_width[:, 0]) - min(lat_max_width[:, 0])
    lat_cc_x_limit = max(n_slice_2D[:, 0]) - 0.35 * slice_width

    med_max_width = n_slice_2D[n_slice_2D[:, 0] <
                               (min(n_slice_2D[:, 0]) + 0.4 * slice_width)]
    med_max_height = max(med_max_width[:, 1]) - min(med_max_width[:, 1])
    med_max_width = max(med_max_width[:, 0]) - min(med_max_width[:, 0])
    med_cc_x_limit = min(n_slice_2D[:, 0]) + 0.35 * slice_width

    holder_med = []
    holder_lat = []
    for i in range(1, 40, 1):
        percentage = i / 100

        n_lateral_plane = n_slice_2D[n_slice_2D[:, 0] >
                                     (max(n_slice_2D[:, 0]) - percentage * slice_width)]

        n_medial_plane = n_slice_2D[n_slice_2D[:, 0] <
                                    (min(n_slice_2D[:, 0]) + percentage * slice_width)]

        med_bf_xc, med_bf_yc, med_bf_r, med_sd = circle_fit.least_squares_circle(n_medial_plane)

        lat_bf_xc, lat_bf_yc, lat_bf_r, lat_sd = circle_fit.least_squares_circle(n_lateral_plane)

        # print(i,'% medial:', med_sd, ' lateral: ', lat_sd)

        if med_bf_r < (med_max_height / 2 + 1) and med_bf_xc < med_cc_x_limit:
            # if med_sd < 0.5:
            holder_med.append([med_bf_xc, med_bf_yc, z_value, med_bf_r, med_sd, i])

        if lat_bf_r < (lat_max_height / 2 + 1) and lat_bf_xc > lat_cc_x_limit:
            # if lat_sd < 0.5:
            holder_lat.append([lat_bf_xc, lat_bf_yc, z_value, lat_bf_r, lat_sd, i])

    holder_med = np.array(holder_med)
    holder_lat = np.array(holder_lat)

    if len(holder_med) > 1:
        holder_med = np.column_stack((holder_med, np.sqrt(holder_med[:, 0] ** 2 + holder_med[:, 1] ** 2)))
        med_shortlist = holder_med

        for i in [6, 1, 0, 6, 1, 0]:
            try:
                temp_a, temp_b = np.histogram(med_shortlist[:, i], bins='doane')
                ind = np.where(temp_a == max(temp_a))[0]
                med_shortlist = med_shortlist[med_shortlist[:, i] <= temp_b[ind + 1], :]
                med_shortlist = med_shortlist[med_shortlist[:, i] >= temp_b[ind], :]
                if len(med_shortlist) < 3:
                    break
            except:
                pass

        med_shortlist = med_shortlist[med_shortlist[:, 4] == min(med_shortlist[:, 4])].flatten()
        output_med = [med_shortlist[[0, 1, 2]], med_shortlist[3], med_shortlist[[4, 5]]]
    elif len(holder_med) == 1:
        print(holder_med)
        output_med = [holder_med[0, [0, 1, 2]], holder_med[0, 3], holder_med[0, [4, 5]]]
    else:
        output_med = None

    if len(holder_lat) > 1:
        holder_lat = np.column_stack((holder_lat, np.sqrt(holder_lat[:, 0] ** 2 + holder_lat[:, 1] ** 2)))
        lat_shortlist = holder_lat
        for i in [6, 1, 0, 3, 6, 1, 0, 3]:
            try:
                temp_a, temp_b = np.histogram(lat_shortlist[:, i], bins='doane')
                ind = np.where(temp_a == max(temp_a))[0]
                lat_shortlist = lat_shortlist[lat_shortlist[:, i] <= temp_b[ind + 1], :]
                lat_shortlist = lat_shortlist[lat_shortlist[:, i] >= temp_b[ind], :]
                if len(lat_shortlist) < 3:
                    break
            except:
                pass

        lat_shortlist = lat_shortlist[lat_shortlist[:, 4] == min(lat_shortlist[:, 4])].flatten()
        output_lat = [lat_shortlist[[0, 1, 2]], lat_shortlist[3], lat_shortlist[[4, 5]]]

    elif len(holder_lat) == 1:
        output_lat = [holder_lat[0, [0, 1, 2]], holder_lat[0, 3], holder_lat[0, [3, 4, 5]]]
    else:
        output_lat = None

    return output_med, output_lat


def plateau_parameters_new():
    """temporarily stored here for now"""
    s_plateau = SimpleVtkExtraction(s_in).get_plane_view(ccc, cutPlane=[0, 0, 1])
    s_ML = SimpleVtkExtraction(s_plateau).get_plane_view(ccc, cutPlane=[0, 1, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_ML.GetBounds()
    pl_ML = [xmax - xmin, [xmin, 0, 0], [xmax, 0, 0]]
    s_AP = SimpleVtkExtraction(s_plateau).get_plane_view(ccc, cutPlane=[1, 0, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_AP.GetBounds()
    pl_AP = [ymax - ymin, [0, ymin, 0], [0, ymax, 0]]
    s_MAP = SimpleVtkExtraction(s_in).get_plane_view(mcc, cutPlane=[1, 0, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_MAP.GetBounds()
    pl_MAP = [ymax - ymin, [0, ymin, 0], [0, ymax, 0]]
    s_LAP = SimpleVtkExtraction(s_in).get_plane_view(lcc, cutPlane=[1, 0, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_LAP.GetBounds()
    pl_LAP = [ymax - ymin, [0, ymin, 0], [0, ymax, 0]]
    return pl_ML, pl_AP, pl_MAP, pl_LAP