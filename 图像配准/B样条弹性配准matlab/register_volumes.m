function [Ireg,O_trans,Spacing,Tx,Ty,Tz] = register_volumes(Imoving,Istatic,type)
% This function register_volumes is the most easy way to register two
% 3D images both rigid and nonrigidly.
%
% Features:
% - It can be used with 3D images from different type of scans or modalities.
% - It uses both a rigid transform and a nonrigid b-spline grid transform.
% - It uses grid refinement
% - It can be used with images of different sizes.
% - The function will automaticaly detect if the images can be registered
% with the sum of squared pixel distance (SSD), or when mutual information 
% must be used as image similarity measure.
%
% [Ireg,Grid,Spacing,Tx,Ty,Tz] = register_volumes(Imoving,Istatic,Type);
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
%         3D image in the same way: I=bspline_transform(Grid,I,Spacing); 
%   Spacing: The uniform b-spline knot spacing 
%   Tx, Ty, Tz : The transformation fields of the non-rigid transform in 
%       x,y and z direction seen from the  static image to the moving image.
%
% Example,
%
% % Get the volume data
% [Imoving,Istatic]=get_example_data;
%
% % Register the images
% Ireg = register_volumes(Imoving,Istatic);
%
% % Show the results
% showcs3(Imoving);
% showcs3(Istatic);
% showcs3(Ireg);
%
% Function is written by D.Kroon University of Twente (October 2008)

% Set verbose 0,1,2 (Display information)
verbose = 2;

% Check for presence of needed functions
if(exist('bspline_transform_3d_double')~=3)
    error('bspline_transform_3d_double mex function not found, compile the c-file');
end
if(exist('bspline_transform_3d_single')~=3)
    error('bspline_transform_3d_single mex function not found, compile the c-file');
end
if(exist('rigid_transform_3d_double')~=3)
    error('rigid_transform_3d_double mex function not found, compile the c-file');
end
if(exist('rigid_transform_3d_single')~=3)
    error('rigid_transform_3d_single mex function not found, compile the c-file');
end
if(exist('mutual_histogram_double')~=3)
    error('mutual_histogram_double mex function not found, compile the c-file');
end
if(exist('mutual_histogram_single')~=3)
    error('mutual_histogram_single mex function not found, compile the c-file');
end

% Store the class of the inputs
Iclass=class(Imoving);

% Convert uint8, uint32 etc. to single.
if(~strcmpi(Iclass,'single')&&~strcmpi(Iclass,'double'))
    range=getrangefromclass(Imoving);
    Imoving=single(Imoving)./range(2);
    range=getrangefromclass(Istatic);
    Istatic=single(Istatic)./range(2);
end

% Resize the moving image to fit the static image
if(sum(size(Istatic)-size(Imoving))~=0)
    Imoving=imresize3d(Imoving,[],size(Istatic),'cubic');
end

% Smooth both images for faster registration
ISmoving=imfilter(Imoving,fspecial3('gaussian',[6 6 6]));
ISstatic=imfilter(Istatic,fspecial3('gaussian',[6 6 6]));

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
    

% Parameter caling of the Translation, Resize and Rotation
scale=[1 0.01 1];

if(verbose>0), disp('Start Rigid registration'); drawnow; end
% Rigid register the smoothed images to get the registration parameters

% Minimizer parameters

% Use struct because expanded optimset is part of the Optimization Toolbox.
if (verbose<2), optim=struct('GradObj','off','Display','off','MaxIter',100,'MaxFunEvals',1000);
else optim=struct('GradObj','off','Display','iter','MaxIter',100,'MaxFunEvals',1000); end
x=fminsd(@(x)rigid_registration_error(x,scale,ISmoving,ISstatic,type),[0 0 0 100 100 100 0 0 0],optim);
 
% Scale the translation, resize and rotation parameters to the real values
x=x.*[scale(1) scale(1) scale(1) scale(2) scale(2) scale(2) scale(3) scale(3) scale(3)];

% Make the affine transformation matrix
M=make_transformation_matrix(x(1:3),x(4:6),x(7:9));

% Do the registration
Imoving=rigid_transform(Imoving,M);


% Non-rigid b-spline grid registration
if(verbose>0), disp('Start non-rigid b-spline grid registration'); drawnow; end

% Calculate max refinements steps
MaxItt=min(floor(log2(size(Imoving)/3)));

% set b-spline grid spacing in x and y direction
Spacing=[2^MaxItt 2^MaxItt 2^MaxItt];

% Remove som refinements steps because a grid with one pixel spacing is not needed
MaxItt=MaxItt-2;

% Make the Initial b-spline registration grid
O_trans=make_init_grid(Spacing,size(Imoving));

% Smooth image for fast registration
ISmoving=imfilter(Imoving,fspecial3('gaussian',[6 6 6]));

% Make the Initial b-spline registration grid
O_trans=make_init_grid(Spacing,size(Imoving));

if (verbose>0), disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow; end

% Smooth both images for faster registration
ISmoving=imfilter(Imoving,fspecial3('gaussian',[6 6 6]));
ISstatic=imfilter(Istatic,fspecial3('gaussian',[6 6 6]));

% Optimizer parameters
if(verbose>0),
    optim=struct('Display','iter','GradObj','on','MaxIter',20,'DiffMinChange',0.1,'DiffMaxChange',1);
else
    optim=struct('Display','off','GradObj','on','MaxIter',20,'DiffMinChange',0.1,'DiffMaxChange',1);
end


% set registration options.
options.type=type;
% Enable forward instead of central gradient incase of error measure is pixel distance
if(strcmpi(type,'sd')), options.centralgrad=false; end

% Reshape O_trans from a matrix to a vector.
sizes=size(O_trans); O_trans=O_trans(:);

% Start the b-spline nonrigid registration optimizer
O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,ISmoving,ISstatic,options),O_trans,optim);

% Reshape O_trans from a vector to a matrix
O_trans=reshape(O_trans,sizes);


for refine_itt=1:MaxItt
    if (verbose>0), disp('Registration Refinement'); drawnow; end

    % No smoothing in last registration step
    if(refine_itt==MaxItt), ISmoving=Imoving; ISstatic=Istatic; end
    
    % Refine the b-spline grid
    O_trans=refine_grid(O_trans); Spacing=Spacing/2;
    if (verbose>0), disp(['Current Grid size : ' num2str(size(O_trans,1)) 'x' num2str(size(O_trans,2)) 'x' num2str(size(O_trans,3)) ]); drawnow; end
    
    % Reshape O_trans from a matrix to a vector.
    sizes=size(O_trans); O_trans=O_trans(:);

    % Start the b-spline nonrigid registration optimizer
    O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,ISmoving,ISstatic,options),O_trans,optim);

    % Reshape O_trans from a vector to a matrix
    O_trans=reshape(O_trans,sizes);
end

% Transform the input image with the found optimal grid.
if ( nargout<4 )
    Ireg=bspline_transform(O_trans,Imoving,Spacing);
else
    [Ireg,Tx,Ty,Tz]=bspline_transform(O_trans,Imoving,Spacing);
end

% Set the class of output to input class
if(strcmpi(Iclass,'uint8')), Ireg=uint8(Ireg*((2^8)-1)); end
if(strcmpi(Iclass,'uint16')), Ireg=uint16(Ireg*((2^16)-1)); end
if(strcmpi(Iclass,'uint32')), Ireg=uint32(Ireg*((2^32)-1)); end
if(strcmpi(Iclass,'int8')), Ireg=int8(Ireg*((2^7)-1)); end
if(strcmpi(Iclass,'int16')), Ireg=int16(Ireg*((2^15)-1)); end
if(strcmpi(Iclass,'int32')), Ireg=int32(Ireg*((2^31)-1)); end