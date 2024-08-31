This repository contains extracts from software developed in 2021-2022 for a research. This is a extract only to demonstrate my capabilities, they will not run without the other dependencies which have been omitted on purpose 

# Background
This project first came into idea when my supervisors went to some conferences for orthopaedic surgeries. There were famous orthopaedic surgeons from many countries there, and at some point there was a heated debate on how to approach Total Knee Replacement Surgeries. Surgeons from America would say they need to start with a certain section of the femur because they were the biggest part that obstructed view, while surgeons from France would not understand their viewpoint at all. After few hours of the discussion, the participants at the conference realised the leg bones they would deal with during TKR surgeries are different. The participants at this conference took away that there are variation in shapes and sizes of lower leg bones of humans according to ethnicity and population groups, but there haven't been a proper study done to investigate this properly.

This project aimed to investigate and prove that the variation in shapes and sizes of lower leg bones between populations does exist. If significant difference did exist, than it would advocate another hypothesis that success rate of TKR can be increased by appropriately reflecting the ethnic background of the patient and using different surgical approaches to tailor to their bone shapes.

To acheive the object, codes in this repository was developed to automate the process of measuring different parameters such as dimensions and angles that are used by academics to study the leg bones.
The below code was used to process 464 femur and tibia bones and measured 19 different parameters on femurs and tibia. 

The resultsing measurements were used to perform stastical analysis. Statistical  indicated difference between population groups does exist.

# About the Code
The codes in this repository is intended to only showcase my skills and experience in development, and will not run properly as it has other owned dependeices that was omitted on purpose.

The program written for this research used framework developed in-house that utilised Python implementation of the VTK (Visualisation Toolkit) library that was used to read the 3D model of the bones to allow manipulation. 

The loaded bone was represented in 2D array format that represents each coordinates of the of the 3D surface mesh.

Numpy library was used to work with this corrdinates to identify bone landmarks and perform measurements.

## _01_execution_related.py 
contains the execution script. It handles the process of loading the STL file and converting it to VTK object for handling, identifying landmarks on bones and calling functions from _02_femur_functions and _02_tibia_functions to measure the bone parameters and output the result to a .csv file

## _02_femur_functions.py 
contains functions that are intended to be used for STL files of femur bones only. Functions here will receive VTK object that are segments of the femur and identify landmarks on the femur to perform measurements 

## _02_tibia_functions.py 
contains functions to measure tibia parameters. 

## _03_misc_functions.py 
contains functions used by femur and tibia measurement functtions. These misc functions uses mathematical methods such as circle fitting that are needed to find parameters on the bones. Uses scipy, matplotlib  
