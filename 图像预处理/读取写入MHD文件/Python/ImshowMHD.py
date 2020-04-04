#-*- coding:utf-8 -*-
# Author: Lajipeng
# Date: 2019/8/19
# Fuction:imshow .mhd 
import cv2
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
# cv2 is incompatible with 64-bits opencv-python
# path = r"E:\Research\Liver 20181228 reg\#2\051_YE JIN SHENG_2661199_20180529_028Y_M\YE JIN SHENG_2661199_20180529_028Y_M_001\result.0.mhd"
path = r"C:\Users\10446\Desktop\result.mhd"
image =sitk.ReadImage(path)
# image = sitk.GetArrayFromImage(image)
print(image)
