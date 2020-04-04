# -*- coding:utf-8 -*-
# Author: Lajipeng
# Data: 2019/10/18
# Fuction: Make .bat for controlling elastix which will used for registration(affine)
import os
import sys
import SimpleITK as sitk
import numpy as np
import pydicom

# elastix -f "E:/Research/Liver 20181228 mhd/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000/YE JIN SHENG_2661199_20180529_028Y_M_000.mhd" -m "E:/Research/Liver 20181228 mhd/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_001/YE JIN SHENG_2661199_20180529_028Y_M_001.mhd" -out "E:/Research/Liver 20181228 reg/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000" -p "E:/Research/Registration2/elastix/example/exampleinput/parameters_BSpline.txt"


filename = 'elastix_affine1_.bat'

Text1 = "elastix -f"

Text2 = " -p \"G:/Research/Registration2/elastix/example/exampleinput/parameters_Affine.txt\""


SaveRawDicom = "G:/Research/Liver 20181228 mhd/Liver_cancer2"  # where to save mhd and raw files

SaveRegDicom = "G:/Research/Liver 20181228 reg/Liver_cancer2_affine"  # where to save mhd and raw files

Series_ID = []
with open(filename, 'a') as file_object:
    file_object.write('G:\ncd G:/research/Registration2/elastix/elastix-4.9.0-win64' + '\n')
    F1 = os.listdir(SaveRawDicom)
    for i in F1:
        F2 = os.listdir(SaveRawDicom + '/' + i)
        F1_i = i
        fix_list = []
        moving_list = []
        for j in F2:
            F1_j = SaveRawDicom + '/' + F1_i + '/' + j + '/' + j.replace(' ', '_') + '.mhd'
            F3_j = SaveRegDicom + '/' + i + '/'
            fix_list.append((F1_j,j))
            # sys.exit(0)
        count = 0
        fix_mhd_path = fix_list[0][0]
        fix_list.pop(0)
        fix_list.pop(-1)
        for moving_mhd_path in fix_list:
            file_object.write(Text1 + ' ' + "\"" + str(fix_mhd_path) + "\"" + ' -m ' + "\"" + str(moving_mhd_path[0]) + "\"" + ' -out ' + "\"" + F3_j + str(moving_mhd_path[1]) + "\"" + Text2 + '\n')





