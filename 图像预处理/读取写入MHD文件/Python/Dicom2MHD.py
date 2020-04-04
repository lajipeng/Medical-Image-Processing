#-*- coding:utf-8 -*-
# Author: Lajipeng
# Data: 2019/8/12
# Fuction: Transform dcm files to mhd


import os

import pydicom

import numpy

import SimpleITK as sitk

import sys

from adjust_window import window_adjust

 
#%%
# Path Statement

PathDicom = "E:/Research/Liver 20181228/#3"		# where to save dicom files 

SaveRawDicom = "E:/Research/Liver 20181228 mhd/#3"     # where to save mhd and raw files

lstFilesDCM = []
F1 = os.listdir(PathDicom)
count = 0
for i in F1:
    count = count + 1
    if count <= 2:
        continue    
    F2 = os.listdir(PathDicom + '/' + i)
    F1_i = i
    s = []
    for j in F2:
        F1_j = PathDicom + '/' + F1_i + '/' + j 
        F2_j = SaveRawDicom + '/' + F1_i + '/' + j + '/' + j.replace( ' ' , '_' ) + '.mhd'
        
        # load the path of saving dicom files 

        for dirName, subdirList, fileList in os.walk(F1_j):

                for filename in fileList:

                        if ".dcm" in filename.lower():  # judge dicom files

                                #print(filename)

                                lstFilesDCM.append(F1_j + '/' + filename)  # put it into list
                                
 
#%%
        # 1：The first dcm will be the reference and all dcms have the same dimensions

        RefDs = pydicom.read_file(lstFilesDCM[0])  # read the first dicom image
        

         
#%%
        # 2：Get the dimensions of dicom image

        ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM)) # ConstPixelDimsis a tuple

#%%      

        # 3：Get the Spacing along x and y and the SliceThick along z

        ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

         
#%%
        # 4：Get the Origin of dicom image

        Origin = RefDs.ImagePositionPatient

        # Creat a numpy 3D-array and set data_type

        ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)  # array is a numpy array

         
#%%
        # 5: Read all the dicom image in the list and put the pixel_array in the 3D-array

        
        for filenameDCM in lstFilesDCM:
                #print(filenameDCM)
                ds = pydicom.read_file(filenameDCM)
                #print(ds)
                pixel_array = window_adjust(ds)
                #print(pixel_array)
                ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = pixel_array
        lstFilesDCM.clear()
        

         
#%%
        # 6：Transpose the numpy array

        ArrayDicom = numpy.transpose(ArrayDicom, (2, 0, 1))
#%%
        # 7: Transform the numpy array into mhd and raw through SimpleITK

        sitk_img = sitk.GetImageFromArray(ArrayDicom, isVector=False)

        sitk_img.SetSpacing(ConstPixelSpacing)

        sitk_img.SetOrigin(Origin)

        sitk.WriteImage(sitk_img, F2_j)

        print("%s successfully transformed!" % F1_j)

        #sys.exit(0)
