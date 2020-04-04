#-*- coding:utf-8 -*-
# Author: Lajipeng
# Data: 2019/8/12
# Fuction: Make .bat for controlling elastix which will used for registration
import os
import sys
import SimpleITK as sitk
import numpy as np
import pydicom
# elastix -f "E:/Research/Liver 20181228 mhd/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000/YE JIN SHENG_2661199_20180529_028Y_M_000.mhd" -m "E:/Research/Liver 20181228 mhd/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_001/YE JIN SHENG_2661199_20180529_028Y_M_001.mhd" -out "E:/Research/Liver 20181228 reg/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000" -p "E:/Research/Registration2/elastix/example/exampleinput/parameters_BSpline.txt"
# transformix -in "E:/Research/Liver 20181228 mhd/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_Seg_001/YE JIN SHENG_2661199_20180529_028Y_M_Seg_001.mhd" -out "E:/Research/Liver 20181228 reg/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000" -tp "E:/Research/Liver 20181228 reg/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000/TransformParameters.0.txt"

filename = 'elastix.bat'

Text1 = "elastix -f"

Text2 = " -p \"H:/Research/Registration2/elastix/example/exampleinput/parameters_BSpline.txt\""

Text3 = "transformix -in"

Text4 = "/TransformParameters.0.txt"

PathDicom = "E:/Research/Liver 20181228/#2"		# where to save dicom files 

SaveRawDicom = "E:/Research/Liver 20181228 mhd/Liver_cancer"     # where to save mhd and raw files

SaveRegDicom = "E:/Research/Liver 20181228 reg/Liver_cancer2"     # where to save mhd and raw files


Series_ID = []
with open (filename,'a') as file_object :
    file_object.write('E:\ncd E:/research/Registration2/elastix/elastix-4.9.0-win64' + '\n')
    F1 = os.listdir(SaveRawDicom)
    for i in F1:
        F2 = os.listdir(SaveRawDicom + '/' + i)
        F1_i = i
        fix_list = []
        for j in F2:
            F1_j = SaveRawDicom + '/' + F1_i + '/' + j + '/' + j.replace( ' ' , '_' ) + '.mhd' 
            F2_j = PathDicom + '/' + i + '/' + j + '/' + '00000001.dcm'
            F3_j = SaveRegDicom + '/' + i + '/'
            RefDs = pydicom.read_file(F2_j)
            fix_list.append((F1_j,RefDs.SeriesNumber,j))
            #sys.exit(0)
        flag = PathDicom + '/' + i
        count = 0
        for fix_mhd in fix_list:
            if fix_mhd[1] == fix_list[-1][1]:
                if fix_mhd == fix_list[-1]:
                    print("Error, no counterpart/n Please input an integer after checking the file:%s\n" % flag)
                    index = input("input")
                    index = int(index)
                    fix_mhd_path = fix_list[index][0]
                    moving_mhd_path = fix_list.pop(index)
                    seg_mhd_path = fix_list.pop(-1)    
                else:
                    moving_mhd_path = fix_mhd[0]
                    fix_list.pop(count)
                    seg_mhd_path = fix_list.pop(-1) 
                    break 
            count = count + 1
        for fix_mhd in fix_list:
            file_object.write(Text1 + ' ' +  "\"" + str(fix_mhd[0]) + "\"" + ' -m ' + "\"" + str(moving_mhd_path) + "\"" + ' -out ' + "\"" + F3_j  + str(fix_mhd[2]) + "\"" + Text2 + '\n' )
            file_object.write(Text3 + ' ' +  "\"" + str(seg_mhd_path[0]) + "\"" + ' -out ' + "\"" + F3_j  + str(fix_mhd[2]) + "\"" + ' -tp ' + "\"" + F3_j + str(fix_mhd[2]) + Text4 + "\"" + '\n' )
            
            
            
       


