# Import Libraries
import numpy as np
import tensorflow as tf
import os
import time
os.environ['CUDA_VISIBLE_DEVICES'] = '5,6'
class L0Smooth_Tf():
    '''
    image shape:[N,W,1]
    image value:[0,1]
    '''    
    def __init__(self,MTF):
        self.kappa = 2.0
        self._lambda= 0.08
        self.MTF = MTF
        self.windowx = tf.reshape(tf.constant([[-1.0,1.0]],dtype=tf.float32),[1,2,1,1,1])
        self.windowy = tf.reshape(tf.constant([[-1.0,1.0]],dtype=tf.float32),[2,1,1,1,1])
        self.windowz = tf.reshape(tf.constant([[-1.0,1.0]],dtype=tf.float32),[1,1,2,1,1])
        self.windowdx = tf.reshape(tf.constant([[1.0,-1.0]],dtype=tf.float32),[1,2,1,1,1])
        self.windowdy = tf.reshape(tf.constant([[1.0,-1.0]],dtype=tf.float32),[2,1,1,1,1])
        self.windowdz = tf.reshape(tf.constant([[1.0,-1.0]],dtype=tf.float32),[1,1,2,1,1])


    def Smooth(self,S):
        # Compute image OTF
        N, M, D = [128,128,128]

        # Compute F(I)
        FI = tf.signal.fft3d(tf.cast(tf.complex(S, tf.zeros([N, M, D],dtype=tf.float32)),dtype = tf.complex64))
        
        # Iteration settings
        beta_max = 1e5
        beta = 2 * self._lambda
        iteration = 0
        # Iterate until desired convergence in similarity
        while beta < beta_max:
            # start_time = time.time()
            ### Step 1: estimate (h, v) subproblem
            # subproblem 1

            # compute dxSp
            padx = tf.reshape(S[:,0,:],[128,1,128])
            h = tf.reshape(tf.concat([S,padx],axis=1),[1,128,129,128,1])
            h = tf.reshape(tf.nn.conv3d(h,self.windowx,[1,1,1,1,1],padding='VALID'),[128,128,128,1])

            # compute dySp
            pady = tf.reshape(S[0,:,:],[1,128,128])
            v = tf.reshape(tf.concat([S,pady],axis=0),[1,129,128,128,1])
            v = tf.reshape(tf.nn.conv3d(v,self.windowy,[1,1,1,1,1],padding='VALID'),[128,128,128,1])
            
            # compute dzSp
            padz = tf.reshape(S[:,:,0],[128,128,1])
            k = tf.reshape(tf.concat([S,padz],axis=2),[1,128,128,129,1])
            k = tf.reshape(tf.nn.conv3d(k,self.windowz,[1,1,1,1,1],padding='VALID'),[128,128,128,1])

            # compute minimum energy E = dxSp^2 + dySp^2 <= _lambda/beta
            t = 3*(tf.pow(h, 2) + tf.pow(v, 2) + tf.pow(k, 2)) >= self._lambda / beta
            t = tf.cast(t, dtype = tf.float32)

            # compute piecewise solution for hp, vp
            h = tf.multiply(h,t)
            v = tf.multiply(v,t)
            k = tf.multiply(k,t)

            ### Step 2: estimate S subproblem
            # subproblem 2 

            # compute dxhp + dyvp
            padx = tf.reshape(h[:,127,:,:],[128,1,128,1])
            dxhp = tf.reshape(tf.concat([padx,h],axis=1),[1,128,129,128,1])
            dxhp = tf.reshape(tf.nn.conv3d(dxhp,self.windowdx,[1,1,1,1,1],padding='VALID'),[128,128,128])
            pady = tf.reshape(v[127,:,:,:],[1,128,128,1])
            dyvp = tf.reshape(tf.concat([pady,v],axis=0),[1,129,128,128,1])
            dyvp = tf.reshape(tf.nn.conv3d(dyvp,self.windowdy,[1,1,1,1,1],padding='VALID'),[128,128,128])
            padz = tf.reshape(k[:,:,127,:],[128,128,1,1])
            dzkp = tf.reshape(tf.concat([padz,k],axis=2),[1,128,128,129,1])
            dzkp = tf.reshape(tf.nn.conv3d(dzkp,self.windowdz,[1,1,1,1,1],padding='VALID'),[128,128,128])
            normin = dxhp + dyvp + dzkp
            normin = tf.cast(tf.complex(normin, tf.zeros(normin.shape,dtype=tf.float32)),dtype=tf.complex64)

            FS = tf.signal.fft3d(normin)

            # solve for S + 1 in Fourier domain
            denorm = tf.constant(1,dtype=tf.complex64) + beta * self.MTF
            FS = (FI + beta * FS) / denorm
            
            # inverse FFT to compute S + 1
            temp = tf.signal.ifft3d(FS)
            S = tf.math.real(temp)
            # update beta for next iteration
            beta *= self.kappa
            iteration += 1
            # final_time = time.time()
            # print("Total Time: %f (s)" % (final_time - start_time))
        # Rescale image
        S = tf.tile(S[tf.newaxis,:,:,:,tf.newaxis],[1,1,1,1,1])
        return S*255