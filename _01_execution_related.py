import glob
import os
import numpy as np
import pandas as pd
from PIL import Image

"""libraries"""
# owned
from _02_femur_functions import *
from _02_tibia_functions import *
from _03_misc_functions import *
from msk_align import SimpleAlignSVDFemur
from msk_math import *
from msk_vtk_fun import *
from msk_vtk import SimpleCoordSystem, SimpleVtkHead, SimpleVtkBody, SimpleVtkTail, \
    ComplexVtkExtraction, SimpleVtkExtraction, \
    SimpleVtkDefaultShape

from vtk.util.colors import cold_grey



def visualise(input_vtk):  # input vtk in the form of [part1,part2,part3]
    """
    add the vtk objects in a list format
    function automatically launches vtk window to visualise the objects
    closing the vtk window terminates the script
    """
    renderer = SimpleCoordSystem().get_renderer()
    renderer = SimpleVtkBody(input_vtk, renderer=renderer, opacity=0.2).get_renderer()
    SimpleVtkTail(renderer=renderer, Position=[0, -23000 / (16 / 9), 0], FocalPoint=[0, 0, 0],
                  Viewup=[0, 0, 1], ParallelProjection=False, Viewangle=(1)).get_interactor()


def screenshots(vtk_list, mode=None, FocalPoint=[0, 0, 0], reverse = False,
                autoscale = True, opacity = 0.5, suf ='', f_in_s ='',aspect_ratio = 16/9):
    """
    takes screenshot of the vtk objects
    automatically measures the size of the object to get the best zoom possible
    function designed for on 1920 x 1080 monitor, using different size may give undesired results
    """
    if cp_single_image_test_mode:
        global f_in

    if f_in_s != '':
        f_in = f_in_s

    try:
        f_in = os.path.join(dir, filename)
    except:
        pass



    if suf != '':
        suf = str(suf)
        len_suf = len(suf)
        zeros = 5-len_suf
        zeros = zeros *'0'
        suf = '_' + zeros + str(suf)+'.stl'
        f_in = f_in.replace('.stl',suf)


    z_diff = 610

    y_diff = 87
    x_diff = 80
    z_midpoint = FocalPoint[2]

    if autoscale:
        if len(vtk_list)>1:
            s_all = vtk_list[0]
            for i in vtk_list:
                s_all = add_append_poly_filter(s_all,i)
        else:
            s_all = vtk_list[0]
        xmin, xmax, ymin, ymax, zmin, zmax = s_all.GetBounds()
        z_diff = zmax-zmin
        if z_diff ==0:
            z_diff = 1
        y_diff = ymax - ymin
        if y_diff ==0:
            y_diff = 1
        x_diff = xmax - xmin
        if x_diff == 0:
            x_diff =1
        z_midpoint = (zmax + zmin) / 2
        y_midpoint = (ymax + ymin) / 2
        x_midpoint = (xmax + xmin) / 2


        FocalPoint = [x_midpoint,y_midpoint,z_midpoint]

    #y_pos = -z_diff/math.tan(math.pi*1/180) / (16 / 9)
    rotate_image = True
    default_sagittal_viewup = [0, 1, 0]
    default_frontal_viewup = [-1, 0, 0]

    if z_diff/x_diff < aspect_ratio:
        cor_pos = -x_diff/math.tan(math.pi*1/180)
        rotate_image = False
        default_frontal_viewup = [0,0,1]

    else:
        cor_pos = -z_diff/math.tan(math.pi*1/180) / aspect_ratio

    if z_diff/y_diff <aspect_ratio:
        sag_pos = -y_diff/math.tan(math.pi*1/180)
        rotate_image = False
        default_sagittal_viewup = [0, 0, 1]

    else:
        sag_pos = -z_diff / math.tan(math.pi * 1 / 180) / aspect_ratio

    axial_z_pos = (y_diff/math.tan(math.pi*1/180))+z_midpoint
    #axial_z_pos = y_diff

    renderer = SimpleCoordSystem().get_renderer()
    renderer = SimpleVtkBody(vtk_list, renderer=renderer, opacity=opacity).get_renderer()

    frontal_pos_default = np.array([0, cor_pos, 0]) + np.array(FocalPoint)
    saggital_pos_default = np.array([sag_pos, 0, 0]) + np.array(FocalPoint)
    axial_pos_default = np.array([0, 0, axial_z_pos]) + np.array(FocalPoint)




    if reverse:
        frontal_pos_default[1] = -frontal_pos_default[1]
        saggital_pos_default[0] = -saggital_pos_default[0]
        axial_pos_default[2] = -axial_pos_default[2]
        f_in = f_in.replace('.stl','_zrv.stl')

    if mode.lower() == 'frontal':
        SimpleVtkTail(renderer=renderer, FocalPoint=FocalPoint, Viewangle=(1),
                      Position=frontal_pos_default, Viewup=default_frontal_viewup,
                      imageSave=f_in.replace('.stl', '_frontal.png', )
                      )
        if rotate_image:
            picture = Image.open(f_in.replace('.stl', '_frontal.png'))

            picture.rotate(90, expand=True).save(f_in.replace('.stl', '_frontal.png'))

    if mode.lower() == 'saggital':
        SimpleVtkTail(renderer=renderer, FocalPoint=FocalPoint, Viewangle=(1),
                      Position=saggital_pos_default, Viewup=default_sagittal_viewup,
                      imageSave=f_in.replace('.stl', '_sagittal.png', )
                      )
        if rotate_image:
            picture = Image.open(f_in.replace('.stl', '_sagittal.png'))

            picture.rotate(90, expand=True).save(f_in.replace('.stl', '_sagittal.png'))

    if mode.lower() == 'axial':
        SimpleVtkTail(renderer=renderer, FocalPoint=FocalPoint, Position=axial_pos_default, Viewup=[0, 1, 0],
                      viewangle=(1), imageSave=f_in.replace('.stl', '_axial.png', ))
    if mode.lower() == 'all':

        if z_diff / x_diff > 16 / 9:
            rotate_image = True

        SimpleVtkTail(renderer=renderer, FocalPoint=FocalPoint, Viewangle=(1), Position=frontal_pos_default,
                      Viewup=default_frontal_viewup,
                      imageSave=f_in.replace('.stl', '_frontal.png', )
                      )
        if rotate_image:

            picture = Image.open(f_in.replace('.stl', '_frontal.png'))

            picture.rotate(90, expand=True).save(f_in.replace('.stl', '_frontal.png'))

        if z_diff / y_diff > 16 / 9:
            rotate_image = True
        else:
            rotate_image = False

        SimpleVtkTail(renderer=renderer, FocalPoint=FocalPoint, Position=saggital_pos_default, Viewup=default_sagittal_viewup
                      , Viewangle=(1), imageSave=f_in.replace('.stl', '_sagittal.png', )
                      )
        if rotate_image:

            picture = Image.open(f_in.replace('.stl', '_sagittal.png'))

            picture.rotate(90, expand=True).save(f_in.replace('.stl', '_sagittal.png'))

        SimpleVtkTail(renderer=renderer, FocalPoint=FocalPoint, Position=axial_pos_default, Viewup=[0, 1, 0],
                      viewangle=(1), imageSave=f_in.replace('.stl', '_axial.png', ))

    if suf != '':

        f_in = f_in.replace(suf,'.stl')

    return


def calculate_and_save_tibia(f_in):
    height, prox_tibia_points, dist_tibia_points, diameter_20, diameter_50, diameter_80, \
    pl_mediolateral_points, pl_ap_points, pl_map_points, pl_lap_points, s_in = get_all_parameters_tibia(
        f_in)

    results = [height[0], prox_tibia_points[2], dist_tibia_points[2], diameter_20[2],
               diameter_50[2], diameter_80[2], pl_mediolateral_points[2], pl_ap_points[2],
               pl_map_points[2], pl_lap_points[2]]

    """below is the screenshot function """
    if cp_take_screenshots:
        renderer = SimpleCoordSystem().get_renderer()

        # create shapes for each parameter
        s_height = SimpleVtkDefaultShape().get_cylinder_by_points(1, [-10, 0, height[1][2]],
                                                                  [-10, 0, height[2][2]])
        s_prox_tibia = SimpleVtkDefaultShape().get_cylinder_by_points(1, prox_tibia_points[0], prox_tibia_points[1])
        s_dist_tibia = SimpleVtkDefaultShape().get_cylinder_by_points(1, dist_tibia_points[0], dist_tibia_points[1])
        s_diameter_20 = SimpleVtkDefaultShape().get_cylinder_by_points(1, diameter_20[0], diameter_20[1])
        s_diameter_50 = SimpleVtkDefaultShape().get_cylinder_by_points(1, diameter_50[0], diameter_50[1])
        s_diameter_80 = SimpleVtkDefaultShape().get_cylinder_by_points(1, diameter_80[0], diameter_80[1])
        s_pl_mediolateral = SimpleVtkDefaultShape().get_cylinder_by_points(0.3, pl_mediolateral_points[0],
                                                                           pl_mediolateral_points[1])
        s_pl_ap = SimpleVtkDefaultShape().get_cylinder_by_points(0.3, pl_ap_points[0], pl_ap_points[1])
        s_pl_map = SimpleVtkDefaultShape().get_cylinder_by_points(0.3, pl_map_points[0], pl_map_points[1])
        s_pl_lap = SimpleVtkDefaultShape().get_cylinder_by_points(0.3, pl_lap_points[0], pl_lap_points[1])

        s_anatomic_axis = SimpleVtkDefaultShape().get_cylinder_by_points(1, prox_tibia_points[3],
                                                                         dist_tibia_points[3])

        s_out = [s_in, s_height, s_prox_tibia, s_dist_tibia, s_diameter_20, s_anatomic_axis,
                 s_diameter_50, s_diameter_80, s_pl_mediolateral, s_pl_ap, s_pl_map, s_pl_lap]

        screenshots(s_out, mode='all')

    return results


def get_all_parameters_tibia(f_in):
    """loading in the tibia and checking whether its a right tibia,
    flips the bone if it is right"""
    s_in = SimpleAlignSVDFemur(f_in, isRight=False).get_aligned_femur()

    (xmax, xmin, ymax, ymin, zmin, zmax) = s_in.GetBounds()

    pos_x_bounds = [xmin, xmax, ymin, ymax, zmin, zmin + 3]
    pos_x_shape_out = SimpleVtkExtraction(s_in).get_cube(pos_x_bounds, inside_out=True)

    _, xval, _, _, _, _ = pos_x_shape_out.GetBounds()

    if xval > 0:

        print('Right Tibia detected')

        s_in = SimpleAlignSVDFemur(f_in, isRight=True).get_aligned_femur()

        if cp_output_csv:
            global filename
            os.rename(f_in, dir + 'zzz' + filename)
            f_in = os.path.join(dir, filename.replace('.stl', '_sub.stl'))
            write(s_in, f_in)
            # right femur mirrored and printed to stl
            filename = filename.replace('.stl', '_sub.stl')

    """ align to TAA"""
    taa, taa_coords, s_shaft = taa_svd(s_in)
    s_taa = SimpleVtkDefaultShape().get_cylinder_by_points(
        0.5, taa_coords[0], taa_coords[1])
    taa_midpoint = (taa_coords[0] + taa_coords[1]) / 2
    print(taa_midpoint)
    print(taa)

    xy_translate_coord = -(taa_coords[0] + taa_coords[1]) / 2

    xy_translate_coord = list(xy_translate_coord)

    s2t = add_transform_by_translation(s_in, xy_translate_coord)

    if taa[2] < 0:
        taa = -taa
    rotate_angle = cal_angle_by_vectors(taa, [0, 0, 1])
    rotate_angle = [rotate_angle[3], rotate_angle[1], rotate_angle[2]]

    if taa[0] > 0:
        rotate_angle[1] = -rotate_angle[1]
    if taa[1] < 0:
        rotate_angle[2] = -rotate_angle[2]
    s2 = add_transform_by_rotation(s2t, rotation=rotate_angle)

    """finished aligning to TAA"""

    points = extract_points_from_shape(s2)

    height = get_height(points)
    prox_ml_poin, _, _, _, _, _ = mediolateral_dimensions(s2)

    """ initial tibia alignment to temporary tibia plateau axis
    temporary tibia plateau axis = furthest points on slice 20mm from tip of tibia"""
    plane_20mm = SimpleVtkExtraction(s2).get_plane_view([0, 0, height[1][2] - 20], [0, 0, 1])
    furthest_point_axis = furthest_points_in_slice(extract_points_from_shape(plane_20mm))

    furthest_point_axis = furthest_point_axis[1] - furthest_point_axis[0]
    if furthest_point_axis[0] < 0:
        furthest_point_axis = -furthest_point_axis
    rotate_angle = cal_angle_by_vectors(furthest_point_axis, [1, 0, 0])
    rotate_angle = [rotate_angle[3], rotate_angle[2], rotate_angle[1]]
    if furthest_point_axis[1] > 0:
        rotate_angle = [-rotate_angle[0], 0, 0]
    s2 = add_transform_by_rotation(s2, rotation=rotate_angle)

    """find the plateau axis and rotate tibia to the axis"""
    xy_axis = tibia_plateau_axis(s2, height, screenshot=False)

    tp_axis = xy_axis[0]
    mag_tp_axis = np.array(tp_axis)
    mag_tp_axis = np.sqrt(np.sum(np.square(mag_tp_axis)))
    if tp_axis[0] < 0:
        tp_axis = -tp_axis
    xy_translate_coord = list(-xy_axis[1])
    cc_z = -xy_translate_coord[2]
    xy_translate_coord[2] = 0
    s2 = add_transform_by_translation(s2, xy_translate_coord)
    rotate_angle = cal_angle_by_vectors(tp_axis, [1, 0, 0])
    rotate_angle = [rotate_angle[3], rotate_angle[1], rotate_angle[2]]
    if tp_axis[1] > 0:
        rotate_angle = [-rotate_angle[0], 0, 0]

    s_in = add_transform_by_rotation(s2, rotation=rotate_angle)
    points = extract_points_from_shape(s_in)
    mcc = np.array([-0.5 * mag_tp_axis, 0, cc_z])
    lcc = np.array([0.5 * mag_tp_axis, 0, cc_z])
    ccc = np.array([0, 0, cc_z])
    """

    End of alignment bit

    """
    height = get_height(points)  # height returns 3 components, index = 0 is height

    prox_ML_midpoint, dist_ML_midpoint, prox_lateral, prox_medial, dist_lateral, dist_medial = mediolateral_dimensions(
        s_in)  # returns points for all
    prox_tibia_points = [prox_lateral, prox_medial, np.sqrt(np.sum((prox_lateral - prox_medial) ** 2, axis=0)),
                         prox_ML_midpoint]
    dist_tibia_points = [dist_lateral, dist_medial, np.sqrt(np.sum((dist_lateral - dist_medial) ** 2, axis=0)),
                         dist_ML_midpoint]

    diameter_20, diameter_50, diameter_80 = tibia_shaft_diameters(s_in, points, height)

    s_plateau = SimpleVtkExtraction(s_in).get_plane_view(ccc, cutPlane=[0, 0, 1])
    s_ML = SimpleVtkExtraction(s_plateau).get_plane_view(ccc, cutPlane=[0, 1, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_ML.GetBounds()
    pl_ML = [[xmin, 0, cc_z], [xmax, 0, cc_z], xmax - xmin]

    s_AP = SimpleVtkExtraction(s_plateau).get_plane_view(ccc, cutPlane=[1, 0, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_AP.GetBounds()
    pl_AP = [[0, ymin, cc_z], [0, ymax, cc_z], ymax - ymin]

    s_MAP = SimpleVtkExtraction(s_plateau).get_plane_view(mcc, cutPlane=[1, 0, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_MAP.GetBounds()
    pl_MAP = [[xmax, ymin, cc_z], [xmax, ymax, cc_z], ymax - ymin]

    s_LAP = SimpleVtkExtraction(s_plateau).get_plane_view(lcc, cutPlane=[1, 0, 0])
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_LAP.GetBounds()
    pl_LAP = [[xmax, ymin, cc_z], [xmax, ymax, cc_z], ymax - ymin]

    return height, prox_tibia_points, dist_tibia_points, diameter_20, diameter_50, diameter_80, \
           pl_ML, pl_AP, pl_MAP, pl_LAP, s_in


def calculate_and_save_femur(f_in):
    """section different parts of the femur and return them as vtk shapes"""
    s_in, s_fh, s_les_trochanter, s_proxima, s_shaft, s_condyle, s_shaft_top, s_shaft_bot = \
        read_and_section_femur(f_in)
    (xmin, xmax, ymin, ymax, zmin, zmax) = s_in.GetBounds()
    output_result = []

    l_height = SimpleVtkDefaultShape().get_cylinder_by_points(1, [0, 0, zmax], [0, 0, zmin])

    l_anatomic, v_f_anatomic = cal_faa(s_shaft)
    angle_faa_fma = cal_angle_by_vectors([0, 0, 1], v_f_anatomic)
    output_result = output_result + list(angle_faa_fma)

    l_femural_shaft_top, v_femur_shaft_top = cal_fst(s_shaft_top)
    l_femural_shaft_bot, v_femur_shaft_bot = cal_fsb(s_shaft_bot)
    angle_bow = cal_angle_by_vectors(v_femur_shaft_top, v_femur_shaft_bot)
    output_result = output_result + list(angle_bow)

    l_diaphyseal_condyle, v_disphyseal_condyle = cal_diaphyseal_condyle(s_condyle)
    angle_joint_line = cal_angle_by_vectors(v_f_anatomic, v_disphyseal_condyle)
    output_result = output_result + list(angle_joint_line)

    l_femural_neck, v_femur_neck = cal_fna(s_proxima)
    angle_version_TEA = cal_angle_by_vectors([1, 0, 0], v_femur_neck)
    output_result = output_result + list(angle_version_TEA)

    l_posterior_condyle, v_posterior_condyle = cal_posterior_condyle(s_condyle)
    angle_version_pca = cal_angle_by_vectors(v_femur_neck, v_posterior_condyle)
    output_result = output_result + list(angle_version_pca)

    angle_femoral_neck = cal_angle_by_vectors(v_femur_neck, v_f_anatomic)
    output_result = output_result + list(angle_femoral_neck)

    l_mechanical, s_fh_fit, radius, fhc = cal_fhc(s_fh)
    output_result.append(radius)
    l_fh_radius = SimpleVtkDefaultShape().get_cylinder_by_points(1, fhc, fhc + [-radius, 0, 0])

    meidaleteral_langth = cal_meidaleteral_langth(s_condyle)
    output_result.append(meidaleteral_langth)

    femur_langth = cal_femur_length(s_in)
    output_result.append(femur_langth)

    s_tea = SimpleVtkExtraction(s_in).get_plane_view([0, 0, 0], [0, 0, 1])
    s_tea = SimpleVtkExtraction(s_tea).get_plane_view([0, 0, 0], [0, 1, 0])
    xmin, xmax, _, _, _, _ = s_tea.GetBounds()
    l_tea = SimpleVtkDefaultShape().get_cylinder_by_points(2, [xmin, 0, 0], [xmax, 0, 0])
    if cp_take_screenshots:
        # s_out = [s_in,s_fh_fit]
        # screenshots(s_out,mode='all',opacity=0.3,aspect_ratio=2000/2000)
        # os._exit(0)
        s_out = add_append_poly_filter(l_anatomic, l_mechanical, l_femural_neck, l_femural_shaft_top,
                                       l_femural_shaft_bot, l_diaphyseal_condyle, l_posterior_condyle,
                                       s_fh_fit)
        s_all = [s_in, s_fh, s_les_trochanter, s_proxima, s_shaft, s_condyle, s_shaft_top, s_shaft_bot, s_out]
        z_mid = np.mean([zmax, zmin])
        x_mid = np.mean([xmax, xmin])
        y_mid = np.mean([ymax, ymin])

        s_out = add_append_poly_filter(l_anatomic, l_mechanical, l_femural_neck, l_femural_shaft_top,
                                       l_femural_shaft_bot, l_diaphyseal_condyle, l_posterior_condyle,
                                       )

        screenshots(s_all, mode='all', autoscale=True, opacity=0.8)
        screenshots(s_all, mode='all', reverse=True, autoscale=True, opacity=0.8)

        pass

    return output_result


def read_and_section_femur(f_in):
    # --------------------------------------------------------------------------------------
    # read the shape, center of mass and aignment

    s_in = SimpleAlignSVDFemur(f_in, isRight=False, plot_check=False).get_aligned_femur()

    # --------------------------------------------------------------------------------------
    # extract femoral head

    (xmin, xmax, ymin, ymax, zmin, zmax) = s_in.GetBounds()

    """
    detect if femur is right femur by checking the top most point of z axis is 
    higher than the zmax of the shape
    """

    pos_x_bounds = [xmax - 30, xmax, ymin, ymax, zmax - 40, zmax]
    pos_x_shape_out = SimpleVtkExtraction(s_in).get_cube(pos_x_bounds, inside_out=True)

    neg_x_bounds = [xmin, xmin + 30, ymin, ymax, zmax - 40, zmax]
    neg_x_shape_out = SimpleVtkExtraction(s_in).get_cube(neg_x_bounds, inside_out=True)

    _, _, _, _, _, pos_z = pos_x_shape_out.GetBounds()
    _, _, _, _, _, neg_z = neg_x_shape_out.GetBounds()

    """Detect if femur head is on the left or right, if on right flips the and save the mirrored image for future use"""
    if pos_z > neg_z:
        print('right femur detected')
        s_in = SimpleAlignSVDFemur(f_in, isRight=True).get_aligned_femur()

        (xmin, xmax, ymin, ymax, zmin, zmax) = s_in.GetBounds()

        if cp_output_csv:
            global filename
            os.rename(f_in, dir + 'zzz' + filename)
            f_in = os.path.join(dir, filename.replace('.stl', '_sub.stl'))
            write(s_in, f_in)
            filename = filename.replace('.stl', '_sub.stl')

    s_fh = ComplexVtkExtraction(s_in).get_femoral_head()

    """
    identify lesser trochanter
    """
    points_list = extract_points_from_shape(SimpleVtkExtraction(s_in).get_cube_by_range([0, 100]))
    points_list = points_list[
        points_list[:, 0] > max(points_list[:, 0]) - (max(points_list[:, 0]) - min(points_list[:, 0])) * 0.2]

    """greater trochanter identified and 7cm below is taken to be searched for lesser trochanter"""
    z_g_trochanter = points_list[points_list[:, 1] == max(points_list[:, 1]), 2]
    z_offset_trochanter = z_g_trochanter - 70  # lower point of less trochanter

    s_les_trochanter = SimpleVtkExtraction(s_in).get_cube_by_range(
        clip_range=[abs(zmin) + z_offset_trochanter, zmax - z_g_trochanter],
        inside_out=True)  # cut off the trochanter + 20mm

    points_list = extract_points_from_shape(s_les_trochanter)

    """
    Loop from greater trochanter to bottom to find the lower trochanter
    lesser trochanter found by comparing X and Y coordinate to section 5mm above the point in question
    (lesser trochanter protudes so assumption is if x or y value is greater 5mm above it its the trochanter)
    variable to hold the z coordinate of the lesser trochanter"""
    z_coord_less_trochanter = points_list[points_list[:, 0] ==
                                          min(points_list[points_list[:, 2] <
                                                          max(points_list[:, 2]) - 10, 0]), 2]
    points_temp = points_list
    while True:

        points_5_up = points_list[points_list[:, 2] > z_coord_less_trochanter]
        points_5_up = points_5_up[points_5_up[:, 2] < z_coord_less_trochanter + 5]

        points_temp = points_temp[points_temp[:, 2] <= z_coord_less_trochanter, :]
        xmin_points_list = points_temp[points_temp[:, 0] == min(points_temp[:, 0]), 0]
        ymin_points_list = points_temp[points_temp[:, 0] == max(points_temp[:, 0]), 1]
        xmin_point_up = points_5_up[points_5_up[:, 0] == min(points_5_up[:, 0]), 0]
        ymin_point_up = points_5_up[points_5_up[:, 0] == max(points_5_up[:, 0]), 1]

        if xmin_point_up > xmin_points_list:
            break

        if ymin_points_list > ymin_point_up:
            if abs(xmin_points_list - xmin_point_up) < 0.5:
                break
        """scans the next 5mm"""
        points_temp = points_temp[points_temp[:, 2] <= z_coord_less_trochanter - 5, :]

        z_coord_less_trochanter = points_temp[points_temp[:, 0] == min(points_temp[:, 0]), 2]

    # --------------------------------------------------------------------------------------
    # get proximal part of the femur
    s_proxima = SimpleVtkExtraction(s_in).get_cube_by_range(
        clip_range=[abs(zmin) + z_coord_less_trochanter, 0], inside_out=True)
    s_shaft = SimpleVtkExtraction(s_in).get_cube_by_range(
        clip_range=[abs(zmin), zmax - z_coord_less_trochanter], inside_out=True)
    s_condyle = SimpleVtkExtraction(s_in).get_cube_by_range(
        clip_range=[abs(zmin), 0], inside_out=False)

    (_, _, _, _, zmin_shaft, zmax_shaft) = s_shaft.GetBounds()
    shaft_length = zmax_shaft - zmin_shaft
    chop_off_perc = 0.1  # chop off the outer 10% of the shape to avoid the outliers
    s_shaft_bot = SimpleVtkExtraction(s_shaft).get_cube_by_range(
        clip_range=[chop_off_perc * shaft_length, 0.5 * shaft_length], inside_out=True)
    s_shaft_top = SimpleVtkExtraction(s_shaft).get_cube_by_range(
        clip_range=[0.5 * shaft_length, chop_off_perc * shaft_length], inside_out=True)

    return s_in, s_fh, s_les_trochanter, s_proxima, s_shaft, s_condyle, s_shaft_top, s_shaft_bot


def function_to_run():
    """insert test function here"""
    return


if __name__ == '__main__':

    """
    control panel
    """
    # running mode
    cp_single_image_test_mode = False
    cp_multiple_image_test_mode = False
    cp_take_screenshots = True
    cp_output_csv = True
    cp_bone_mode = 'tibia'
    cp_raise_error = True

    # directory specification
    dir = 'C:/Users/Sky/Desktop/FYP - SSM/Shuqiao feb ver/Measure_Auto/00_stls'

    # test file if test mode
    f_in = './00_stls/surface_scan/f_002_lol_07.stl'

    # f_in = 'C:/Users/Sky/Desktop/FYP - SSM/STL from SQ harddrive/NMDID/From ruonan group/bin/TL/m_121343_TL.stl'

    files = os.listdir(dir)
    print(files)

    """
    end of control panel
    """

    errors = []
    data_dic = {}
    if cp_single_image_test_mode:
        if f_in == '':
            print('specify f_in directory and filename')
            exit()
        """
        refer the code to test here"""
        function_to_run(f_in)
        exit()
    for filename in files:
        need_mirroring = False
        if '.stl' in filename.lower():
            if 'zzz' in filename.lower():
                continue
            try:
                f_in = os.path.join(dir, filename)
                print(f_in)
                # if '_R.stl'.lower() in filename.lower():
                #     need_mirroring = True
                #
                # if 'FR.stl'.lower() in filename.lower():
                #     need_mirroring = True
                #
                """
                add code from here
                """
                if cp_multiple_image_test_mode:
                    function_to_run(f_in)
                    # temp_func(f_in)
                    continue

                if cp_bone_mode.lower() == 'tibia':
                    data_dic[filename] = calculate_and_save_tibia(f_in)
                    pass

                if cp_bone_mode.lower() == 'femur':
                    print('femur mode')
                    data_dic[filename] = calculate_and_save_femur(f_in)
                    print('no error yet')

                """
                stop here
                """
            except:
                if cp_raise_error:
                    raise Exception("dfsdsdsd")
                data_dic[filename] = "error"
                print('error written to dic')
                errors.append(filename)
                pass
    print('Files with errors: ', errors)
    if cp_output_csv:

        if cp_bone_mode.lower() == 'tibia':
            output_column = ['Height', 'Prox ML Width', 'Dist ML Width',
                             'Diameter at 20%', 'Diameter at 50%', 'Diameter at 80%',
                             'Plateau ML Width', 'Plateau AP Length',
                             'Plateau MAP Length', 'Plateau LAP Length']
            df_data = pd.DataFrame.from_dict(
                data_dic, orient='index', columns=output_column)

            df_data.to_csv('output.csv', mode='a', sep=',')

        if cp_bone_mode.lower() == 'femur':
            output_column = ['FAA - FMA Angle', 'cor', 'sag', 'trans',
                             'Bow Angle', 'cor', 'sag', 'trans',
                             'Joint Line - FAA Angle', 'cor', 'sag', 'trans',
                             'Femoral neck - Mediolateral Angle', 'cor', 'sag', 'trans',
                             'Femoral neck - posterior line angle', 'cor', 'sag', 'trans',
                             'Femoral neck - Anatomical Axis', 'cor', 'sag', 'trans',
                             'Femoral Head Radius', 'Mediolateral Length', 'Femur Length']
            df_data = pd.DataFrame.from_dict(
                data_dic, orient='index', columns=output_column)

            df_data.to_csv('output.csv', mode='a', sep=',')
