import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv
class tsmooth():
#    tsmooth - Structure Extraction from Texture via Relative Total Variation
#    S = tsmooth(I, lambda, sigma, maxIter) extracts structure S from
#    structure+texture input I, with smoothness weight lambda, scale
#    parameter sigma and iteration number maxIter.

#    Paras:
#    @I         : Input UINT8 image, both grayscale and color images are acceptable.
#    @lambda    : Parameter controlling the degree of smooth.
#                 Range (0, 0.05], 0.01 by default.
#    @sigma     : Parameter specifying the maximum size of texture elements.
#                 Range (0, 6], 3 by defalut.f
#    @sharpness : Parameter controlling the sharpness of the final results,
#                 which corresponds to \epsilon_s in the paper [1]. The smaller the value, the sharper the result.
#                 Range (1e-3, 0.03], 0.02 by defalut.
#    @maxIter   : Number of itearations, 4 by default.
#    Example
#    ==========
#    I  = imread('Bishapur_zan.jpg');
#    S  = tsmooth(I); % Default Parameters (lambda = 0.01, sigma = 3, sharpness = 0.02, maxIter = 4)
#    figure, imshow(I), figure, imshow(S);

#    ==========
#    The Code is created based on the method described in the following paper
#    [1] "Structure Extraction from Texture via Relative Total Variation", Li Xu, Qiong Yan, Yang Xia, Jiaya Jia, ACM Transactions on Graphics,
#    (SIGGRAPH Asia 2012), 2012.
#    The code and the algorithm are for non-comercial use only.

#    Author: Wang Peng (wangp16@fudan.edu.com)
#    Date  : 03/09/2020
#    Version : 1.0
#    Copyright 2020, Fudan University.
    def __init__(self, I, lbd, sigma, sharpness, maxIter):
        self.lbd = 0.01
        self.sigma = 3.0
        self.sharpness = 0.02
        self.maxIter = 4
        self.I = np.asarray(I/255.0, dtype=float)
        if lbd is not None:
            self.lbd = lbd
        if sigma is not None:
            self.sigma = sigma
        if sharpness is not None:
            self.sharpness = sharpness
        if maxIter is not None:
            self.maxIter = maxIter
    def RTV(self):
        x = self.I
        sigma_iter = self.sigma
        lamb = self.lbd/2.0
        dec=2.0
        for iter in range(0,self.maxIter):
            [wx, wy] = self.computeTextureWeights(x, sigma_iter)
            x = self.solveLinearEquation(self.I, wx, wy, lamb)
            sigma_iter = sigma_iter/dec
            if sigma_iter < 0.5:
                sigma_iter = 0.5
        S = x
        return S

    def padarray(self, f, axis, size):
        if axis == -1:
            pad = np.zeros([size[0], size[1], 1])
            ret = np.concatenate((f, pad),axis = axis)
        else:
            pad = np.zeros([size[0],1,size[2]])
            ret = np.concatenate((f, pad),axis = axis)
        return ret

    def lpfilter(self, FImg, sigma):
        FBImg = FImg
        for ic in range(0,FBImg.shape[2]):
            FBImg[:,:,ic] = self.conv2_sep(FImg[:,:,ic], sigma)
        return FBImg

    def conv2_sep(self, im, sigma):
        ksize = np.bitwise_or(int(np.round(5*sigma)),1)
        g = self.Gaussian_filter(ksize, sigma)
        g = np.reshape(g,(1,ksize))
        ret = self.conv2(im,g)
        ret = self.conv2(ret,np.transpose(g))
        return ret

    def conv2(self, A, B):
        ma, na = A.shape
        mb, nb = B.shape
        row_pad = mb -1
        column_pad = nb - 1
        mc = ma + 2*row_pad - mb + 1
        nc = na + 2*column_pad - nb + 1
        A = np.pad(A, ((row_pad,row_pad),(column_pad,column_pad)), 'constant')
        B = np.flip(B)
        C = np.zeros([mc,nc])
        for i in range(0,mc):
            for j in range(0,nc):
                C[i,j] = np.sum(np.multiply(A[i:i+mb,j:j+nb],B))
        ret = C[int(np.ceil((mc-ma)/2)):int(np.ceil((mc+ma)/2)),int(np.ceil((nc-na)/2)):int(np.ceil((nc+na)/2))]
        return ret

    def Gaussian_filter(self, size, sigma):
        eps = 2.2204e-16
        siz = (size-1)/2
        std = sigma
        x = np.linspace(-siz, siz, size)
        arg = -np.square(x)/(2*std*std)
        h = np.exp(arg)
        h = np.maximum(h, eps*np.max(h))
        
        sumh = np.sum(h)
        if sumh != 0:
            h  = h/sumh
        return h

        
    def computeTextureWeights(self, x, sigma):
        fx = np.diff(x, n = 1, axis = 0)
        fx = np.pad(fx,((0,1),(0,0),(0,0)),'constant')
        fy = np.diff(x, n = 1, axis = 1)
        # fy = padarray(fy, axis = -2, size = fy.shape)
        fy = np.pad(fy,((0,0),(0,1),(0,0)),'constant')
        vareps_s = self.sharpness
        vareps = 0.001
        # x1 = np.square(fx)+np.square(fy)
        # x2 = np.sqrt(x1)
        # x3 = np.sum(x2,2)/3
        # x4 = np.max(x3,vareps_s)
        # wto = np.power(x4,-1)
        wto = np.power(np.maximum(np.sum(np.sqrt(np.square(fx)+np.square(fy)),2)/3,vareps_s),-1)
        fbin = self.lpfilter(x, sigma)
        gfy = np.diff(fbin,n=1,axis=0)
        gfy = np.pad(gfy,((0,1),(0,0),(0,0)),'constant')
        gfx = np.diff(fbin,n=1,axis=1)
        gfx = np.pad(gfx, ((0,0),(0,1),(0,0)), 'constant')
        wtbx = np.power(np.maximum(np.sum(np.abs(gfx),2)/3,vareps),-1)
        wtby = np.power(np.maximum(np.sum(np.abs(gfy),2)/3,vareps),-1)
        retx = np.multiply(wtbx,wto)
        rety = np.multiply(wtby,wto)

        retx[:,-1] = 0
        rety[-1,:] = 0
        return [retx, rety]

    def solveLinearEquation(self, IN, wx, wy, lbd):
        #  The code for constructing inhomogenious Laplacian is adapted from
        #  the implementaion of the wlsFilter.
        #  For color images, we enforce wx and wy be same for three channels
        #  and thus the pre-conditionar only need to be computed once.
        (r,c,ch) = IN.shape
        k = r*c
        wx = np.reshape(np.transpose(wx),k)
        dx = -lbd*wx
        wy = np.reshape(np.transpose(wx),k)
        dy = -lbd*wy
        B = np.zeros([2,k])
        B[0,:] = dx
        B[1,:] = dy
        d = np.asarray([-r,-1])
        A = spdiags(B, d, k, k)
        e = dx
        w = np.pad(dx, (r,0), 'constant')
        w = w[0:-r]
        s = dy
        n = np.pad(dy, (1,0), 'constant')
        n = n[0:-1]
        D = 1-(e+w+s+n)
        D = np.reshape(D,[1,k])
        t = np.asarray([0])
        h = spdiags(D, t, k, k)
        A = A + np.transpose(A) + h
        A = csc_matrix(A, dtype=float)
        OUT = np.zeros([r,c,ch],dtype=float)
        for ii in range(0,ch):
            tin = np.reshape(np.transpose(IN[:,:,ii]),[k,1])
            tin = csc_matrix(tin, dtype=float)
            tout = spsolve(A,tin,permc_spec='MMD_AT_PLUS_A')
            temp = np.reshape(tout, [r, c])
            OUT[:,:,ii] = temp
        return OUT
    
# Demo script
# Uncomment each case to see the results
from PIL import Image 
import matplotlib.pyplot as plt
import time
import tensorflow as tf 
t1 = time.time()
I = Image.open('imgs/261.jpg')
I = np.asarray(I)
S = tsmooth(I,0.015,3.0,0.02,4)
result = S.RTV()
t2 = time.time()
print(t2-t1)
# result = Image.fromarray(result)
plt.figure(2)
ax = plt.subplot(1,2,1)
plt.imshow(I)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax = plt.subplot(1,2,2)
plt.imshow(result)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.show()


                    

    