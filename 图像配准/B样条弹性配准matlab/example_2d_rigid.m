% example of 2D affine registration using the lsqnonlin optimizer.

% clean
clear all; close all; clc;

% Read two greyscale images of Lena
I1=im2double(imread('lenag3.png')); 
I2=im2double(imread('lenag2.png'));

% Type of registration error used see registration_error.m
type='d';

% Smooth both images for faster registration
I1s=imfilter(I1,fspecial('gaussian'));
I2s=imfilter(I2,fspecial('gaussian'));

% Parameter caling of the Translation, Resize and Rotation
scale=[1 0.01 1];

[x]=lsqnonlin(@(x)rigid_registration_image(x,scale,I1s,I2s,type),[0 0 100 100 0],[],[],optimset('Display','iter','MaxIter',100));

% Scale the translation, resize and rotation parameters to the real values
x=x.*[scale(1) scale(1) scale(2) scale(2) scale(3)];

% Make the affine transformation matrix
M=make_transformation_matrix(x(1:2),x(3:4),x(5));

Icor=rigid_transform(I1,M);

% Show the registration results
figure,
subplot(1,3,1), imshow(I1);
subplot(1,3,2), imshow(I2);
subplot(1,3,3), imshow(Icor);
