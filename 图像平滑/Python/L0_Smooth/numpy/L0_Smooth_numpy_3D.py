# Import Libraries
import numpy as np
class L0Smooth():
    '''
    image shape:[N,W,1]
    image value:[0,1]
    '''

    def __init__(self, MTF):
        self.kappa = 2.0
        self._lambda = 0.08
        self.MTF = MTF

    def Smooth(self, S):
        # Compute image OTF
        N, M, D = [128, 128, 128]

        # Compute F(I)
        FI = np.complex64(np.zeros((N, M, D, 1)))
        FI[:, :, :, 0] = np.fft.fftn(S[:, :, :, 0])

        # Initialize buffers
        h = np.float32(np.zeros((N, M, D, 1)))
        v = np.float32(np.zeros((N, M, D, 1)))
        k = np.float32(np.zeros((N, M, D, 1)))
        dxhp = np.float32(np.zeros((N, M, D, 1)))
        dyvp = np.float32(np.zeros((N, M, D, 1)))
        dzkp = np.float32(np.zeros((N, M, D, 1)))
        FS = np.complex64(np.zeros((N, M, D, 1)))

        # Iteration settings
        beta_max = 1e5
        beta = 2 * self._lambda
        iteration = 0
        # Iterate until desired convergence in similarity
        while beta < beta_max:
            ### Step 1: estimate (h, v) subproblem
            # subproblem 1

            # compute dxSp
            h[:, 0:M - 1, :, :] = np.diff(S, 1, 1)
            h[:, M - 1:M, :, :] = S[:, 0:1, :, :] - S[:, M - 1:M, :, :]

            # compute dySp
            v[0:N - 1, :, :, :] = np.diff(S, 1, 0)
            v[N - 1:N, :, :, :] = S[0:1, :, :, :] - S[N - 1:N, :, :, :]

            # compute dzSp
            k[:,:,0:D-1,:] = np.diff(S, 1, 2)
            k[:,:,D-1:D,:] = S[:,:,0:1,:] - S[:,:,D-1:D,:]# S[:,:,0:1,:] - S[:,:,D-1:D,:]

            # compute minimum energy E = dxSp^2 + dySp^2 <= _lambda/beta
            t = 3 * (np.power(h, 2) + np.power(v, 2) + np.power(k, 2)) < self._lambda / beta

            # compute piecewise solution for hp, vp
            h[t] = 0
            v[t] = 0
            k[t] = 0

            ### Step 2: estimate S subproblem
            # subproblem 2

            # compute dxhp + dyvp
            dxhp[:, 0:1, :, :] = h[:, M - 1:M, :] - h[:, 0:1, :]
            dxhp[:, 1:M, :, :] = -(np.diff(h, 1, 1))
            dyvp[0:1, :, :, :] = v[N - 1:N, :, :] - v[0:1, :, :]
            dyvp[1:N, :, :, :] = -(np.diff(v, 1, 0))
            dzkp[:,:,0:1,:] = k[:, :, D-1:D] - k[:, :, 0:1]
            dzkp[:,:,1:D,:] = -(np.diff(k, 1, 2))
            normin = dxhp + dyvp + dzkp

            FS[:, :, :, 0] = np.fft.fftn(normin[:, :, :, 0])

            # solve for S + 1 in Fourier domain
            denorm = 1 + beta * self.MTF
            FS[:, :, :, :] = (FI + beta * FS) / denorm
            # print(FS)

            # inverse FFT to compute S + 1
            temp = np.fft.ifftn(FS[:, :, :, 0])
            S[:, :, :, 0] = np.float32(temp.real)
            # update beta for next iteration
            beta *= self.kappa
            iteration += 1
        return S
