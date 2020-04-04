function M=make_transformation_matrix(t,s,r)
% This function make_transformation_matrix.m creates an affine 
% 2D or 3D transformation matrix from translation, resize, and rotation parameters
%
% M=make_transformation_matrix.m(t,s,r)
%
% inputs (2D),
%   t: vector [translateX translateY]
%   s: vector [resizeX resizeY]
%   r: vector [rotate]
% inputs (3D),
%   t: vector [translateX translateY translateZ]
%   s: vector [resizeX resizeY resizeZ]
%   r: vector [rotateX rotateY rotateZ]
% outputs,
%   M: 2D or 3D affine transformation matrix
%
% example,
%   M=transformation_matrix([0.5 0 0],[1 1 1.2],[0 0 0])
% 
% Function is written by D.Kroon University of Twente (October 2008)

% Calculate affine transformation matrix
if(length(r)==1)
    M=mat_tra_2d(t)*mat_siz_2d(s)*mat_rot_2d(r); 
else
    M=mat_tra_3d(t)*mat_siz_3d(s)*mat_rot_3d(r); 
end


function M=mat_rot_2d(r)
	M=[ cos(r) sin(r) 0;
	   -sin(r) cos(r) 0;
	   0 0 1];
   
function M=mat_siz_2d(s)
	M=[1/s(1) 0 0;
	   0 1/s(2) 0;
	   0 0 1];
   
function M=mat_tra_2d(t)
	M=[1 0 -t(1);
	   0 1 -t(2);
	   0 0 1];


function M=mat_rot_3d(r)
    r=r*(pi/180);
    Rx=[1 0 0 0;
        0 cos(r(1)) -sin(r(1)) 0;
        0 sin(r(1)) cos(r(1)) 0;
        0 0 0 1];

    Ry=[cos(r(2)) 0 sin(r(2)) 0;
        0 1 0 0;
        -sin(r(2)) 0 cos(r(2)) 0;
        0 0 0 1];

    Rz=[cos(r(3)) -sin(r(3)) 0 0;
        sin(r(3)) cos(r(3)) 0 0;
        0 0 1 0;
        0 0 0 1];
    M=Rx*Ry*Rz;

function M=mat_siz_3d(s)
	M=[1/s(1) 0 0 0;
	   0 1/s(2) 0 0;
	   0 0 1/s(3) 0;
	   0 0 0 1];

function M=mat_tra_3d(t)
	M=[1 0 0 -t(1);
	   0 1 0 -t(2);
	   0 0 1 -t(3);
	   0 0 0 1];
