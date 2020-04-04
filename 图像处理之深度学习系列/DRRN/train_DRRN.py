import NCC
from keras.layers import Input, Dense, Conv2D, MaxPooling2D, UpSampling2D
from keras.layers import Activation, add,ReLU
from keras.models import Model
from keras import optimizers
from keras import backend as K
import tensorflow as tf
#%% Model Structure
input_img = Input(shape=(128, 128, 1))  # adapt this if using `channels_first` image data format

residual = input_img
R0 = ReLU()(input_img)
inputs = Conv2D(1, (3, 3), strides=1 , padding='same')(R0)
# inputs = self.input(self.relu(x)) 
out = inputs
for _ in range(9):
    R1 = ReLU()(out)
    conv1 = Conv2D(128, (3, 3), strides=1 , padding='same')(R1)
    R2 = ReLU()(conv1)
    out = Conv2D(128, (3, 3), strides=1 , padding='same')(R2)
    out = add([out, inputs])
    # out = self.conv2(self.relu(self.conv1(self.relu(out))))
    # out = torch.add(out, inputs)
R3 = ReLU()(out)
out = Conv2D(1, (3, 3), strides=1 , padding='same')(R3)
# out = self.output(self.relu(out))
decoded = add([out, residual])
# out = torch.add(out, residual)

autoencoder = Model(input_img, decoded)
adam = optimizers.Adam(lr=0.0002, decay=1e-4)

def RelativeTotalVariation_loss(input_img, decoded):
    # batch_size = img1.shape.as_list()[0]
    loss_gen = NCC.NCC_RTV_loss(input_img, decoded, lbd=None, sigma=None, sharpness=None)
    loss = loss_gen.ComputeRTV_loss()
    return loss

autoencoder.compile(optimizer='adam', loss=RelativeTotalVariation_loss)
autoencoder.summary()

#%% Load Data
import numpy as np
from PIL import Image

(x_train, x_test) = np.load("./obj.npy",allow_pickle=True)
x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_train = np.reshape(x_train, (len(x_train), 128, 128, 1))  # adapt this if using `channels_first` image data format
x_test = np.reshape(x_test, (len(x_test), 128, 128, 1))  # adapt this if using `channels_first` image data format
print (x_train.shape)
print (x_test.shape)

#%% Start trainning and save the model
from keras.callbacks import TensorBoard

autoencoder.fit(x_train, x_train,
                epochs=20,
                batch_size=128,
                shuffle=True,
                validation_data=(x_test, x_test),
                callbacks=[TensorBoard(log_dir='/tmp/autoencoder')])

autoencoder.save('my_model.h5')
#%% Visualize the result
import matplotlib.pyplot as plt
decoded_imgs = autoencoder.predict(x_test[0:8])
n = 8  # how many digits we will display
plt.figure(figsize=(20, 5), dpi=100)
for i in range(n):
    # display original
    ax = plt.subplot(2, n, i + 1)
    plt.imshow(x_test[i].reshape(128, 128))
    plt.gray()
    ax.get_xaxis().set_visible(True)
    ax.get_yaxis().set_visible(False)
    
    # SSIM Encode
    ax.set_title("Encode_Image")
    
    npImg = x_test[i]
    npImg = npImg.reshape((128,128))
    formatted = (npImg * 255 / np.max(npImg)).astype('uint8')
    img = Image.fromarray(formatted)

    # display reconstruction
    ax = plt.subplot(2, n, i + 1 + n)
    plt.imshow(decoded_imgs[i].reshape(128, 128))
    plt.gray()
    ax.get_xaxis().set_visible(True)
    ax.get_yaxis().set_visible(False)

plt.show()