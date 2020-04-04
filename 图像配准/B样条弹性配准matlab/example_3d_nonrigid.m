% Example of 3D non-rigid registration
% using steepest gradient optimizer and grid refinement.

% clean
clear all; close all; clc;

% Get the volume data
[I1,I2]=get_example_data;

% Show the volume data
showcs3(I1);
showcs3(I2);

% Start b-spline grid dimensions (in the image)
Spacing=[26 26 12];

% Type of registration error used see registration_error.m
options.type='sd';
% Fast forward error gradient instead of central gradient.
options.centralgrad='false';

% Make the Initial b-spline registration grid
O_trans=make_init_grid(Spacing,size(I1));

% Convert all values tot type double
I1=double(I1); I2=double(I2); O_trans=double(O_trans); 

% Resize I1 to fit I2
I1=imresize3d(I1,[],size(I2),'linear');

% Smooth both images for faster registration
I1s=imfilter(I1,fspecial3('gaussian',[6 6 6]));
I2s=imfilter(I2,fspecial3('gaussian',[6 6 6]));

% Optimizer parameters
optim=struct('Display','iter','GradObj','on','MaxIter',5,'DiffMinChange',0.1,'DiffMaxChange',1);

% Reshape O_trans from a matrix to a vector.
sizes=size(O_trans); O_trans=O_trans(:);

% Start the b-spline nonrigid registration optimizer
O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,options),O_trans,optim);

% Reshape O_trans from a vector to a matrix
O_trans=reshape(O_trans,sizes);

% Refine the b-spline grid
O_trans=refine_grid(O_trans); Spacing=Spacing/2;

% Reshape O_trans from a matrix to a vector.
sizes=size(O_trans); O_trans=O_trans(:);

% Start the b-spline nonrigid registration optimizer
O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,type),O_trans,optim);

% Reshape O_trans from a vector to a matrix
O_trans=reshape(O_trans,sizes);

% Transform the input image with the found optimal grid.
Icor=bspline_transform(O_trans,I1,Spacing); 

% Make a (transformed) grid image
Igrid=make_grid_image(Spacing,size(I1));
Igrid=bspline_transform(O_trans,Igrid,Spacing); 

% Show the registration results
figure,
showcs3(Icor);
showcs3(Igrid);
