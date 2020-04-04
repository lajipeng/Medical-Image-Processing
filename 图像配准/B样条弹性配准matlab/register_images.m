function [Ireg,O_trans,Spacing,Tx,Ty] = register_images(Imoving,Istatic,type)
% This function register_images is the most easy way to register two
% images both rigid and nonrigidly.
%
% Features:
% - It can be used with images from different type of scans or modalities.
% - It uses both a rigid transform and a nonrigid b-spline grid transform.
% - It uses grid refinement
% - It can be used with images of different sizes.
% - The function will automaticaly detect if the images can be registered
% with the sum of squared pixel distance (SSD), or when mutual information 
% must be used as image similarity measure.
%
% [Ireg,Grid,Spacing,Tx,Ty] = register_images(Imoving,Istatic,Type);
%
% Inputs,
%   Imoving : The image which will be registerd
%   Istatic : The image on which Imoving will be registered
% (Optional)
%   Type : Similarity measure used can be set to sd, mi, d, 
%               gd, gc, cc, pi, ld see image_difference.m.
%
% Outputs,
%   Ireg : The registered moving image
%   Grid: The b-spline controlpoints, can be used to transform another
%           image in the same way: I=bspline_transform(Grid,I,Spacing); 
%   Spacing: The uniform b-spline knot spacing 
%   Tx, Ty : The transformation fields of the non-rigid transform in 
%       x and y direction seen from the  static image to the moving image.
%
% Example,
%   % Read two greyscale images of Lena
%   Imoving=imread('lenag1.png'); 
%   Istatic=imread('lenag3.png');
% 
%   % Register the images
%   [Ireg,O_trans,Spacing,Tx,Ty] = register_images(Imoving,Istatic);
%
%   % Show the registration result
%   figure,
%   subplot(2,2,1), imshow(Imoving); title('moving image');
%   subplot(2,2,2), imshow(Istatic); title('static image');
%   subplot(2,2,3), imshow(Ireg); title('registerd moving image');
%   subplot(2,2,4), imshow(Tx,[]); title('Transf. in x direction');
%
% Function is written by D.Kroon University of Twente (October 2008)

% Set verbose 0,1,2 (Display information)
verbose = 2;

% Start time measurement
if(verbose>0), tic; end

% Check for presence of needed functions
if(exist('bspline_transform_2d_double')~=3)
    error('bspline_transform_2d_double mex function not found, compile the c-file');
end
if(exist('rigid_transform_2d_double')~=3)
    error('rigid_transform_2d_double mex function not found, compile the c-file');
end
if(exist('mutual_histogram_double')~=3)
    error('mutual_histogram_double mex function not found, compile the c-file');
end
% Check for optimization functions
use_fminunc=exist('fminunc')>0;

% Store the class of the inputs
Iclass=class(Imoving);

% Convert the inputs to double
Imoving=im2double(Imoving);
Istatic=im2double(Istatic);

% Resize the moving image to fit the static image
if(sum(size(Istatic)-size(Imoving))~=0)
    Imoving = imresize(Imoving,size(Istatic),'bicubic');
end

% Make smooth images for histogram and fast rigid registration
ISmoving=imfilter(Imoving,fspecial('gaussian',[10 10],2.5));
ISstatic=imfilter(Istatic,fspecial('gaussian',[10 10],2.5));

% Detect if the mutual information or pixel distance can be used as 
% similarity measure. By comparing the histograms.
if(exist('type','var')==0)
    Hmoving= hist(ISmoving(:),60)./numel(Imoving);
    Hstatic = hist(ISstatic(:),60)./numel(Istatic);
    if(sum(log(abs(Hmoving-Hstatic)+1))>0.5), 
        type='mi'; 
        if(verbose>0), disp('Using Mutual information'); drawnow; end
    else
        type='sd';
        if(verbose>0), disp('Using Pixel Distance'); drawnow; end
    end
end
% Register the moving image rigid to the static image

% Parameter caling of the Translation, Resize and Rotation
scale=[1 0.01 1];

if(verbose>0), disp('Start Rigid registration'); drawnow; end
% Rigid register the smoothed images to get the registration parameters

% Minimizer parameters
if(use_fminunc)
    if (verbose<2), optim=optimset('GradObj','off','LargeScale','off','Display','off','MaxIter',100,'MaxFunEvals',1000);
    else optim=optimset('GradObj','off','LargeScale','off','Display','iter','MaxIter',100,'MaxFunEvals',1000); end
    x=fminunc(@(x)rigid_registration_error(x,scale,ISmoving,ISstatic,type),[0 0 100 100 0],optim);
else
    % Use struct because expanded optimset is part of the Optimization Toolbox.
    if (verbose<2), optim=struct('GradObj','off','Display','off','MaxIter',100,'MaxFunEvals',1000);
    else optim=struct('GradObj','off','Display','iter','MaxIter',100,'MaxFunEvals',1000); end
    x=fminsd(@(x)rigid_registration_error(x,scale,ISmoving,ISstatic,type),[0 0 100 100 0],optim);
end
 
% Scale the translation, resize and rotation parameters to the real values
x=x.*[scale(1) scale(1) scale(2) scale(2) scale(3)];

% Make the affine transformation matrix
M=make_transformation_matrix(x(1:2),x(3:4),x(5));

% Do the registration
Imoving=rigid_transform(Imoving,M);

% Non-rigid b-spline grid registration
if(verbose>0), disp('Start non-rigid b-spline grid registration'); drawnow; end

% Calculate max refinements steps
MaxItt=min(floor(log2(size(Imoving)/3)));

% set b-spline grid spacing in x and y direction
Spacing=[2^MaxItt 2^MaxItt];

% Remove som refinements steps because a grid with one pixel spacing is not needed
MaxItt=MaxItt-2;

% Make the Initial b-spline registration grid
O_trans=make_init_grid(Spacing,size(Imoving));
if (verbose>0), disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2))]); drawnow; end

% Smooth image for fast registration
ISmoving=imfilter(Imoving,fspecial('gaussian',[10 10],2.5));

% set registration options.
options.type=type;
% Enable forward instead of central gradient incase of error measure is pixel distance
if(strcmpi(type,'sd')), options.centralgrad=false; end

% Reshape O_trans from a matrix to a vector.
sizes=size(O_trans); O_trans=O_trans(:);

% fminunc hessian uses too much memmory with large number of unknowns
if(numel(O_trans)>1024), use_fminunc=false; end
if(use_fminunc)
    % Optimizer parameters
    if (verbose<2), optim=optimset('GradObj','on','LargeScale','off','Display','off','MaxIter',100,'DiffMinChange',0.1,'DiffMaxChange',1,'MaxFunEvals',1000);
    else optim=optimset('GradObj','on','LargeScale','off','Display','iter','MaxIter',100,'DiffMinChange',0.1,'DiffMaxChange',1,'MaxFunEvals',1000); end

    % Start the b-spline nonrigid registration optimizer
    O_trans = fminunc(@(x)bspline_registration_gradient(x,sizes,Spacing,ISmoving,ISstatic,options),O_trans,optim);
else
    % Use struct because expanded optimset is part of the Optimization Toolbox.
    if (verbose<2), optim=struct('GradObj','on','Display','off','MaxIter',100,'DiffMinChange',0.1,'DiffMaxChange',1,'MaxFunEvals',1000);
    else optim=struct('GradObj','on','Display','iter','MaxIter',100,'DiffMinChange',0.1,'DiffMaxChange',1,'MaxFunEvals',1000); end
    
    % Start the b-spline nonrigid registration optimizer
    O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,ISmoving,ISstatic,options),O_trans,optim);
end
    
% Reshape O_trans from a vector to a matrix
O_trans=reshape(O_trans,sizes);

    
for refine_itt=1:MaxItt
    if (verbose>0), disp('Registration Refinement'); drawnow; end

    % No smoothing in last registration step
    if(refine_itt==MaxItt), ISmoving=Imoving; ISstatic=Istatic; end
    
    % Refine the b-spline grid
    O_trans=refine_grid(O_trans); Spacing=Spacing/2;
    
    if (verbose>0), disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2))]); drawnow; end

    % Reshape O_trans from a matrix to a vector.
    sizes=size(O_trans); O_trans=O_trans(:);

    % fminunc hessian uses too much memmory with large number of unknowns
    if(numel(O_trans)>1024), use_fminunc=false; end
    if(use_fminunc)
        % Start the b-spline nonrigid registration optimizer
        O_trans = fminunc(@(x)bspline_registration_gradient(x,sizes,Spacing,ISmoving,ISstatic,options),O_trans,optim);
    else
        % Start the b-spline nonrigid registration optimizer
        O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,ISmoving,ISstatic,options),O_trans,optim);
    end
   
    % Reshape O_trans from a vector to a matrix
    O_trans=reshape(O_trans,sizes);
end

% Transform the input image with the found optimal grid.
[Ireg,Tx,Ty]=bspline_transform(O_trans,Imoving,Spacing); 

% Set the class of output to input class
if(strcmpi(Iclass,'uint8')), Ireg=im2uint8(Ireg); end
if(strcmpi(Iclass,'single')), Ireg=im2single(Ireg); end
if(strcmpi(Iclass,'int16')), Ireg=im2int16(Ireg); end
if(strcmpi(Iclass,'uint16')), Ireg=im2uint16(Ireg); end

% End time measurement
if(verbose>0), toc, end

