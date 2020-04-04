function I=rigid_registration_image(par,scale,I1,I2,type)
% This function registration_rigid_image, uses affine transfomation of the
% 3D input volume and calculates the registration error after transformation.
%
% I=rigid_registration_image(parameters,scale,I1,I2,type);
%
% input,
%   parameters (in 2D) : vector of length 5 : [translateX translateY 
%           resizeX resizeY rotate]
%   parameters (in 3D) : vector of length 9 : [translateX translateY translateZ
%           resizeX resizeY resizeZ rotateX rotateY rotateZ]
%   scale: Scaling of the input parameters [scaleTranslate scaleResize scaleRotate] 
%   I1: The 2D/3D image which is affine transformed
%   I2: The second 2D/3D image which is used to calculate the
%       registration error
%   type: The type of registration error used see image_difference.m
%
% outputs,
%   I: An volume image with the registration error between I1 and I2
%
% example,
%   see example_3d_rigid.m
%
% Function is written by D.Kroon University of Twente (July 2008)

if(ndims(I1)==2)
    par=par.*[scale(1) scale(1) scale(2) scale(2) scale(3)];
    M=make_transformation_matrix(par(1:2),par(3:4),par(5));
    I3=rigid_transform(I1,M);
    [t,I] = image_difference(I3,I2,type);
else
    par=par.*[scale(1) scale(1) scale(1) scale(2) scale(2) scale(2) scale(3) scale(3) scale(3)];
    M=make_transformation_matrix(par(1:3),par(4:6),par(7:9));
    I3=rigid_transform(I1,M);
    [t,I] = image_difference(I3,I2,type);
end

