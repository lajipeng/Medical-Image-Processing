% Example of 3D affine registration using the fminunc optimizer.

% clean
clear all; close all; clc;

% Get the volume data
[I1,I2]=get_example_data;

% Convert all volumes from single to double
I1=double(I1); I2=double(I2);

% First resize volume I1 to match size of volume I2
I1=imresize3d(I1,[],size(I2),'linear');

% Show the volume data
showcs3(I1);
showcs3(I2);

% Type of registration error used see registration_error.m
type='sd';

% Smooth both images for faster registration
I1s=imfilter(I1,fspecial3('gaussian'));
I2s=imfilter(I2,fspecial3('gaussian'));

% Parameter caling of the Translation, Resize and Rotation
scale=[1 0.01 1];

[x]=fminunc(@(x)rigid_registration_error(x,scale,I1s,I2s,type),[0 0 0 100 100 100 0 0 0],optimset('Display','iter','MaxIter',100));
%[x]=lsqnonlin(@(x)rigid_registration_image(x,scale,I1s,I2s,type),[0 0 0 100 100 100 0 0 0],[],[],optimset('Display','iter','MaxIter',100));

% Scale the translation, resize and rotation parameters to the real values
x=x.*[scale(1) scale(1) scale(1) scale(2) scale(2) scale(2) scale(3) scale(3) scale(3)];

% Make the affine transformation matrix
M=make_transformation_matrix(x(1:3),x(4:6),x(7:9));

% Transform the input volume
Icor=rigid_transform(I1,M);

% Show the registration results
figure,
showcs3(Icor);
