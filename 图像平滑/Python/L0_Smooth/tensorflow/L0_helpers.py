import numpy as np
import math

# Compute MTF
def ComputeMTF():
    size_2D = [128,128]
    fx = np.int32([[1, -1]])
    fy = np.int32([[1], [-1]])
    otfFx = psf2otf(fx, size_2D)
    otfFy = psf2otf(fy, size_2D)
    MTF = np.power(np.abs(otfFx), 2) + np.power(np.abs(otfFy), 2)
    MTF = np.tile(MTF[:, :, np.newaxis], (1, 1, 1))
    return MTF

# Compute MTF
def ComputeMTF3D():
    size_3D = [128,128,128]
    fx = np.int32([[[1], [-1]]])
    fy = np.int32([[[1]], [[-1]]])
    fz = np.int32([[[1, -1]]])
    otfFx = psf2otf_3D(fx, size_3D)
    otfFy = psf2otf_3D(fy, size_3D)
    otfFz = psf2otf_3D(fz, size_3D)
    MTF = np.power(np.abs(otfFx), 2) + np.power(np.abs(otfFy), 2) + np.power(np.abs(otfFz), 2)
    # MTF = np.tile(MTF[:, :, :,np.newaxis], (1, 1, 1, 1))
    return MTF

# Convert point-spread function to optical transfer function
def psf2otf(psf, outSize=None):
  # Prepare psf for conversion
  data = prepare_psf(psf, outSize)

  # Compute the OTF
  otf = np.fft.fftn(data)

  return np.complex64(otf)

def psf2otf_3D(psf, outSize=None):
  # Prepare psf for conversion
  data = prepare_psf_3D(psf, outSize)

  # Compute the OTF
  otf = np.fft.fftn(data)

  return np.complex64(otf)

def prepare_psf(psf, outSize=None, dtype=None):
  if not dtype:
    dtype=np.float32

  psf = np.float32(psf)

  # Determine PSF / OTF shapes
  psfSize = np.int32(psf.shape)
  if not outSize:
    outSize = psfSize
  outSize = np.int32(outSize)

  # Pad the PSF to outSize
  new_psf = np.zeros(outSize, dtype=dtype)
  new_psf[:psfSize[0],:psfSize[1]] = psf[:,:]
  psf = new_psf

  # Circularly shift the OTF so that PSF center is at (0,0)
  shift = -(psfSize / 2)
  psf = circshift(psf, shift)

  return psf

def prepare_psf_3D(psf, outSize=None, dtype=None):
  if not dtype:
    dtype=np.float32

  psf = np.float32(psf)

  # Determine PSF / OTF shapes
  psfSize = np.int32(psf.shape)
  if not outSize:
    outSize = psfSize
  outSize = np.int32(outSize)

  # Pad the PSF to outSize
  new_psf = np.zeros(outSize, dtype=dtype)
  new_psf[:psfSize[0],:psfSize[1],:psfSize[2]] = psf[:,:,:]
  psf = new_psf

  # Circularly shift the OTF so that PSF center is at (0,0)
  shift = -(psfSize / 2)
  psf = circshift(psf, shift)

  return psf

# Circularly shift array
def circshift(A, shift):
  for i in range(shift.size):
    A = np.roll(A, int(shift[i]), axis=i)
  return A