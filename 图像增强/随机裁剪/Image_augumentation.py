import os, glob
import numpy as np
from sklearn.model_selection import train_test_split
from PIL import Image
import tensorflow as tf
import matplotlib.pyplot as plt

# Own Your Image Directory
img_dir = ("./Samples/")
img_files = glob.glob(img_dir + "*.jpg")
# Setting Image Propertie
width = 128
height = 128
pixels = width * height * 1 # gray scale 

# Load Image
# AutoEncoder does not have to label data 
x = []
for j in range(0,1):
    for i, f in enumerate(img_files):
        img = Image.open(f)
        img = img.convert("L") # gray sclae
        reshaped_image = tf.cast(img,tf.float32)
        distorted_image = tf.random_crop(reshaped_image,[height,width])
        with tf.Session() as sess:
            distorted_image = sess.run(distorted_image)
        distorted_image = distorted_image/256.0
        label = distorted_image
        label = np.reshape(label,[128,128,1])
        label = np.float32(label)
        x.append([distorted_image,label])
        if i % 10 == 0:
            print(i)
x = np.array(x)
(x_train, x_test) = train_test_split(x, shuffle=False, train_size=0.8, random_state=1)

img_list = (x_train, x_test)
np.save("./Smooth_data.npy", img_list)
print("OK", len(x)) 
