This repository contains codes developed during my final year research project. 
The research project aimed to investigate anatomical variation in human leg bones between different population group.
To acheive the object, codes in this repository was developed to automate the process of measuring different parameters such as dimensions and angles that are used by academics to study the leg bones.
The below code was used to process 464 femur and tibia bones and measured 19 different parameters on femurs and tibia. 
The measurements were used for stastical analysis that indicated difference between population groups does exist.

The codes in this library is intended to be used with VTK library which loads 3D models of leg bones STL file format and create VTK objects to allow manipulations actions such as rotation and translation.

VTK object held the model in the form of point cloud holding the 3D coordinates of each points in 2D array

Numpy library was used to interact with the 2D array and perform mathematical analysis 

# _01_execution_related.py 
contains the execution script. It handles the process of loading the STL file and converting it to VTK object for handling, identifying landmarks on bones and calling functions from _02_femur_functions and _02_tibia_functions to measure the bone parameters and output the result to a .csv file

# _02_femur_functions.py 
contains functions that are intended to be used for STL files of femur bones only. Functions here will receive VTK object that are segments of the femur and identify landmarks on the femur to perform measurements 

# _02_tibia_functions.py 
contains functions to measure tibia parameters. 

# _03_misc_functions.py 
contains functions used by femur and tibia measurement functtions. These misc functions uses mathematical methods such as circle fitting that are needed to find parameters. 
