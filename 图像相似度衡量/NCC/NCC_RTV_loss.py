import tensorflow as tf
import numpy as np

class NCC_RTV_loss():
#    Paras:
#    @I         : Input image, shape = [batchsize,H,W,C]
#    @lambda    : Parameter controlling the degree of smooth.
#                 Range (0, 0.05], 0.01 by default.
#    @sigma     : Parameter specifying the maximum size of texture elements.
#                 Range (0, 6], 3 by defalut.f
#    @sharpness : Parameter controlling the sharpness of the final results,
#                 which corresponds to \epsilon_s in the paper [1]. The smaller the value, the sharper the result.
#                 Range (1e-3, 0.03], 0.02 by defalut.
#    @maxIter   : Number of itearations, 4 by default.
#    ==========
#    The Code is created based on the method described in the following paper
#    [1] "Structure Extraction from Texture via Relative Total Variation", Li Xu, Qiong Yan, Yang Xia, Jiaya Jia, ACM Transactions on Graphics,
#    (SIGGRAPH Asia 2012), 2012.
#    The code and the algorithm are for non-comercial use only.
#    Author: Wang Peng (wangp16@fudan.edu.com)
#    Date  : 03/16/2020
#    Version : 1.0
#    Copyright 2020, Fudan University.

    def __init__(self, I, S, lbd, sigma, sharpness):
        self.lbd = 0.02
        self.sigma = 3.0
        self.sharpness = 0.02
        self.I = tf.transpose(I, [2,0,1,3]) 
        self.S = tf.transpose(S, [2,0,1,3])        
        self.S_b = self.S
        # self.I = I
        # self.S = S
        # self.S_b = tf.transpose(S, [2,0,1,3]) 
        self.channel = 3.0

        if lbd is not None:
            self.lbd = lbd
        if sigma is not None:
            self.sigma = sigma
        if sharpness is not None:
            self.sharpness = sharpness

    def ComputeRTV_loss(self):
        Tv_loss = self.RTV()
        NCC_loss = self.NCC_loss()
        with tf.Session() as sess:
            print("Tv_loss:%f"%sess.run(Tv_loss))
            print("NCC_loss:%f"%(sess.run(NCC_loss)))
        Loss = NCC_loss + self.lbd * Tv_loss
        return Loss
    
    def RTV(self):

        fx = self.tf_diff_axis_0(self.S_b)
        fy = self.tf_diff_axis_1(self.S_b)
        alpha_x2 = tf.reduce_sum(tf.square(fx),axis = 0)/self.channel
        alpha_y2 = tf.reduce_sum(tf.square(fy),axis = 0)/self.channel
        vareps_s = self.sharpness
        vareps = 0.001
        wto = tf.pow(tf.maximum(tf.reduce_sum(tf.sqrt(tf.square(fx)+tf.square(fy)),axis=0)/self.channel,vareps_s),-1)

        fbin = self.lpfilter()
        gfx = self.tf_diff_axis_0(fbin)
        gfy = self.tf_diff_axis_1(fbin)
        wtbx = tf.pow(tf.maximum(tf.reduce_sum(tf.abs(gfx),axis=0)/self.channel,vareps),-1)
        wtby = tf.pow(tf.maximum(tf.reduce_sum(tf.abs(gfy),axis=0)/self.channel,vareps),-1)
        retx = tf.multiply(wtbx,wto)
        rety = tf.multiply(wtby,wto)
        Loss_x = tf.multiply(retx, alpha_x2)
        Loss_y = tf.multiply(rety, alpha_y2)
     
        Tv_loss = tf.reduce_mean(Loss_x + Loss_y)
        return Tv_loss
        
    def lpfilter(self):
        ksize = tf.bitwise.bitwise_or(tf.cast(tf.round(5*self.sigma),'int32'),1)
        g = self.Gaussian_filter(ksize, self.sigma)
        g = tf.reshape(g,[1,ksize,1,1])
        g_transpose = tf.reshape(g,[ksize,1,1,1])
        ret = self.conv2(self.S_b,g)
        ret = self.conv2(ret,g_transpose)
        return ret

    def conv2(self, feature_map, filter):
        '''
        feature_map:shape = [batch_size,h,w,channel]
        filter = [h,w,1,1]
        '''
        ret = tf.nn.conv2d(feature_map, filter, strides=(1,1,1,1), padding='SAME')
        return ret

    def Gaussian_filter(self, size, sigma):
        eps = 2.2204e-16
        siz = (size-1)/2
        std = sigma
        x = tf.linspace(-siz, siz, size)
        arg = -tf.square(x)/(2*std*std)
        h = tf.exp(arg)
        h = tf.maximum(h, eps*tf.reduce_max(h))
        sumh = tf.reduce_sum(h)
        if sumh != 0:
            h  = h/sumh
        return h
    def tf_diff_axis_0(self,a):
        partial_x = a[:,:,1:,:]-a[:,:,:-1,:]
        return tf.pad(partial_x,((0,0),(0,0),(0,1),(0,0)),'constant')

    def tf_diff_axis_1(self,a):
        partial_y = a[:,1:,:,:]-a[:,:-1,:,:]
        return tf.pad(partial_y,((0,0),(0,1),(0,0),(0,0)),'constant')

    def NCC_loss(self):
        sizes = np.prod(self.I.shape.as_list()[1:])
        flatten1 = tf.reshape(self.I, [-1, sizes])
        flatten2 = tf.reshape(self.S, [-1, sizes])
        mean1 = tf.reshape(tf.reduce_mean(flatten1, axis=-1),[-1,1])
        mean2 = tf.reshape(tf.reduce_mean(flatten2, axis=-1),[-1,1])
        var1 = tf.reduce_mean(tf.square(flatten1 - mean1), axis=-1)
        var2 = tf.reduce_mean(tf.square(flatten2 - mean2), axis=-1)
        cov12 = tf.reduce_mean(
            (flatten1 - mean1) * (flatten2 - mean2), axis=-1)
        pearson_r = cov12 / tf.sqrt((var1 + 1e-6) * (var2 + 1e-6))

        ncc_loss = 1 - pearson_r
        ncc_loss = tf.reduce_sum(ncc_loss)
        return ncc_loss

# Demo script
# Uncomment each case to see the results
from PIL import Image 
import matplotlib.pyplot as plt
import time
import tensorflow as tf
t1 = time.time()
I = tf.constant(np.asarray(Image.open('imgs/111.jpg')),shape = [128,128,3,1],dtype=tf.float64)
S = tf.constant(np.asarray(Image.open('imgs/222.jpg')),shape = [128,128,3,1],dtype=tf.float64)
loss_gen1 = NCC_RTV_loss(I/255.0, I/255.0, lbd=None, sigma=None, sharpness=None)
loss_gen2 = NCC_RTV_loss(I/255.0, S/255.0, lbd=None, sigma=None, sharpness=None)
loss1 = loss_gen1.ComputeRTV_loss()
loss2 = loss_gen2.ComputeRTV_loss()
with tf.Session() as sess:
    print("Total loss:",sess.run(loss1))
    print("Total loss:",sess.run(loss2))