function O_trans=make_init_grid(Spacing,sizeI)
% This function creates a uniform 2d or 3D b-spline control grid
% O=make_init_grid(Spacing,sizeI)
% 
%  inputs,
%    Spacing: Spacing of the b-spline grid knot vector
%    sizeI: vector with the sizes of the image which will be transformed
%  
%  outputs,
%    O: Uniform b-splinecontrol grid
%
%  example,
%
%    I1 = im2double(imread('lenag1.png')); 
%    O = make_init_grid([8 8],size(I1));
%
%  Function is written by D.Kroon University of Twente (September 2008)
  
if(length(Spacing)==2)
    % Determine grid spacing
    dx=Spacing(1); dy=Spacing(2);

    % Calculate te grid coordinates (make the grid)
    [X,Y]=ndgrid(-dx:dx:(sizeI(1)+(dx*2)-2),-dy:dy:(sizeI(2)+(dy*2)-2));
    O_trans=ones(size(X,1),size(X,2),2);
    O_trans(:,:,1)=X;
    O_trans(:,:,2)=Y;
else
    % Determine grid spacing
    dx=Spacing(1); dy=Spacing(2); dz=Spacing(3);
    
    % Calculate te grid coordinates (make the grid)
    [X,Y,Z]=ndgrid(-dx:dx:(sizeI(1)+(dx*2)-2),-dy:dy:(sizeI(2)+(dy*2)-2),-dz:dz:(sizeI(3)+(dz*2)-2));
    O_trans=ones(size(X,1),size(X,2),size(X,3),3);
    O_trans(:,:,:,1)=X;
    O_trans(:,:,:,2)=Y;
    O_trans(:,:,:,3)=Z;
end
    
    