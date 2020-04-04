#-- Author：Wang, Peng. Fudan University 
#-- Date：2019/8/13
#-- Fuction：Ajust window_width

import pydicom

import numpy as np

#import sys 


def mat2gray(img,limits):    
        rows = len(img)
        cols = len(img[0])
        for i in range(0,rows):
                for j in range(0,cols):
                        if img[i][j] < 1000:
                                img[i][j] = limits[0]
                                # print(limits[0])
                        elif img[i][j] > limits[1]:
                                img[i][j] = 1064
                                # print(limits[1])
                        else:
                            img[i][j] = 1064
        return img
                                        
def window_adjust(ds):
        # The WindowCenter and WindowWidth are fixed, and there is no need to compute.
        # size = sys.getsizeof(ds.WindowCenter)
        # if size <= 34:
        #         WindowCenter = ds.WindowCenter[0]/ds.RescaleSlope - ds.RescaleIntercept
        #         WindowWidth = ds.WindowWidth[0]/ds.RescaleSlope
        # else:
        #         WindowCenter = ds.WindowCenter/ds.RescaleSlope - ds.RescaleIntercept
        #         WindowWidth = ds.WindowWidth/ds.RescaleSlope
        WindowCenter = 1064.0 
        WindowWidth = 400.0
        limits=[]
        limits.append(WindowCenter - (WindowWidth/2))
        limits.append(WindowCenter + (WindowWidth/2))
        I = mat2gray(ds.pixel_array,limits)
        I = np.uint8((I.astype(float) - 864)/400*255)
        return I

if __name__ == "__main__":
        
        ds = pydicom.read_file('E:/Research/Liver 20181228/#2/051_YE JIN SHENG_2661199_20180529_028Y_M/YE JIN SHENG_2661199_20180529_028Y_M_000/00000001.DCM')
        I = window_adjust(ds)
        print(I)
        
        

                                
                        

