
##QUESTION1:Homography is a matrix which corresponds reference image points to warped image.
##Homography matrix needs both warped image and reference images width and height values to
##compute exact values of the warped image

##QUESTION2: Inverse homography provides us the intencity values. Warped image size And
## reference image size is not same so it is not possible to have one to one relationship.
##We use inverse homography to find intencities of the reference image and fill the Warped
## image with that intencity values.

##QUESTION3: Bilinear interpolation is used for avaraging the points (which are around corresponding)
##point) according to their distance. The pixels with small distances have more effect which is better
##for resolution. In addition this make our image smoother. If we would use closest neighborhood Which
##will take just one pixel. Our output would have worst results.
import cv2
import numpy as np
import math
from numpy.linalg import inv
LAMBDA = 0.5
ROTATION_ANGLE = 30 * math.pi / 180
TILT_AMOUNT = 50 * math.pi / 180

def affineMatrix(rotationAngle, tiltAmount):
    rotationMatrix = np.array([[math.cos(rotationAngle), math.sin(
        rotationAngle)], [-math.sin(rotationAngle), math.cos(rotationAngle)]])
    tilt = np.array([[1/math.cos(tiltAmount), 0.0], [0.0, 1.0]])
    return LAMBDA * (rotationMatrix @ tilt)


def computeSizeofOutput(affineDeformationMatrix, originalShape):
    height, width = originalShape[:2]

    corner1 = affineDeformationMatrix @ np.array([height, 0.0])
    corner2 = affineDeformationMatrix @ np.array([0.0,width])
    corner3 = affineDeformationMatrix @ np.array([height,width])
    corner4 = np.array([0.0, 0.0])
    corners = np.vstack([corner4, corner3, corner2, corner1])
    print(corners)
    xs = []
    ys = []
    ys.append(corner1[0])
    ys.append(corner2[0])
    ys.append(corner3[0])
    ys.append(corner4[0])
    xs.append(corner1[1])
    xs.append(corner2[1])
    xs.append(corner3[1])
    xs.append(corner4[1])

    maxX = max(xs)
    minX = min(xs)
    maxY = max(ys)
    minY = min(ys)
    newWidth = maxX - minX
    newHeight = maxY - minY
    return newHeight, newWidth


def findHomography(affineMatrix, originalShape, newWidth, newHeight):
    first = np.array(
        [[1.0, 0.0, newHeight / 2], [0.0, 1.0, newWidth / 2], [0.0, 0.0, 1.0]])
    second = np.array([[affineMatrix[0][0], affineMatrix[0][1], 0.0], [
                      affineMatrix[1][0], affineMatrix[1][1], 0.0], [0.0, 0.0, 1.0]])
    third = np.array([[1.0, 0.0, -originalShape[0] / 2],
                      [0.0, 1.0, -originalShape[1] / 2], [0.0, 0.0, 1.0]])
    return first @ second @ third


def bilinearInterpolate(img, x, y, height, width):
    minX = int(x)
    minY = int(y)
    if(minX+1 >= width):
        minX = minX - width
    if(minY+1>= height):
        minY = minY -height
    minCenterX = minX + 0.25
    minCenterY = minY + 0.25
    maxCenterX = minX + 0.75
    maxCenterY = minY + 0.75

    intensity2 = ((x - minCenterX) / (maxCenterX - minCenterX)) * img[minY][minX] + (
        (maxCenterX - x) / (maxCenterX - minCenterX)) * img[minY][minX + 1]
    intensity1 = ((x - minCenterX) / (maxCenterX - minCenterX)) * img[minY + 1][minX] + (
        (maxCenterX - x) / (maxCenterX - minCenterX)) * img[minY + 1][minX + 1]
    realIntensity = ((y - minCenterY) / (maxCenterY - minCenterY)) * \
        intensity1 + ((maxCenterY - y) /
                      (maxCenterY - minCenterY)) * intensity2
    return realIntensity


def warpImage(img, invH, newWidth, newHeight):
    dst = np.zeros((int(round(newHeight)), int(round(newWidth))))
    height, width = img.shape[:2]

    i = 0
    while(i < newHeight - 1):
        j = 0
        while(j < newWidth - 1):
            homogPixel = invH @ np.vstack((i, j, 1.0))
            refX = homogPixel[1][0]
            refY = homogPixel[0][0]
            if(refX < 0 or refX > width or refY < 0 or refY > height):
                dst[i][j] = 0
            else:
                dst[i][j] = int(round(bilinearInterpolate(
                    img, refX, refY, height, width)))
            j += 1
        i += 1
    return dst


img = cv2.imread("img1.png", 0)
affine = affineMatrix(ROTATION_ANGLE, TILT_AMOUNT)

newHeight, newWidth = computeSizeofOutput(affine, img.shape)
homography = findHomography(affine, img.shape[:2], newWidth, newHeight)
print(homography)
invH = inv(homography)

dst = warpImage(img, invH, newWidth, newHeight)
cv2.imwrite("result.png", dst)
