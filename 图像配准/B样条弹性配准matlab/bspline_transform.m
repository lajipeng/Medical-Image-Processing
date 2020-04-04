function [I,Tx,Ty,Tz]=bspline_transform(O,I,Spacing)
% Function bspline_transform, is a wrapper of the mex 
% bspline_transform_2d_double and bspline_transform_3d mex functions
%
% I = bspline_transform(O,I,nx,ny,nz)
% I :  Input image.
% O  : Transformation grid 
% Spacing : Are the b-spline grid knot spacings.
%
% Function is written by D.Kroon University of Twente (September 2008)

% init
Tx=0; Ty=0; Tz=0;

% Check if spacing has integer values
if(sum(Spacing-floor(Spacing))>0), error('Spacing must be a integer'); end

if(ndims(I)==2)
    if(~isa(I,'double')), I=im2double(I); end
    if(nargout > 1 )
        [I,Tx,Ty]=bspline_transform_2d_double(O(:,:,1),O(:,:,2),I,Spacing(1),Spacing(2));
    else
        I=bspline_transform_2d_double(O(:,:,1),O(:,:,2),I,Spacing(1),Spacing(2));
    end
else
    if(isa(I,'double'))
        if(nargout > 1 )
            [I,Tx,Ty,Tz]=bspline_transform_3d_double(double(O(:,:,:,1)),double(O(:,:,:,2)),double(O(:,:,:,3)),double(I),double(Spacing(1)),double(Spacing(2)),double(Spacing(3)));
        else
            I=bspline_transform_3d_double(double(O(:,:,:,1)),double(O(:,:,:,2)),double(O(:,:,:,3)),double(I),double(Spacing(1)),double(Spacing(2)),double(Spacing(3)));
        end
    else
        if(nargout > 1 )
             [I,Tx,Ty,Tz]=bspline_transform_3d_single(single(O(:,:,:,1)),single(O(:,:,:,2)),single(O(:,:,:,3)),single(I),single(Spacing(1)),single(Spacing(2)),single(Spacing(3)));
        else
             I=bspline_transform_3d_single(single(O(:,:,:,1)),single(O(:,:,:,2)),single(O(:,:,:,3)),single(I),single(Spacing(1)),single(Spacing(2)),single(Spacing(3)));
        end
   end
end
