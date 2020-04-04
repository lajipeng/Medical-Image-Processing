function [t,I]=image_difference(V,U,type)
% This function [t,I]=image_difference(I1,I2,type) gives the registration error t 
% and error image I between the to images or volumes I1 and I2
% 
% if type,
% 'd'  : differences between I1 and I2
% 'sd' : squared differences
% 'mi' : normalized mutual information
% 'gd' : gradient differences
% 'gc' : gradient correlation
% 'cc' : normalized cros correlation
% 'pi' : pattern intensity
% 'ld' : log absolute difference
% 
%  Example,
%    I1=im2double(imread('lenag1.png')); 
%    I2=im2double(imread('lenag2.png'));
%    [t,I] = registration_error(I1,I2,'sd');
%    disp(t);
%    imshow(I,[])
%
% This function is written by D.Kroon University of Twente (August 2008)

if(exist('type','var')==0), type='p'; end
if(isempty(V)), t=2; I=2; return; end
switch(type)
    case 'd'
        [I,t]=registration_error_differences(V,U);
    case 'sd'
        [I,t]=registration_error_squared_differences(V,U);
    case 'mi'
        [I,t]=registration_error_mutual_info(V,U);
    case 'gd'
        [I,t]=registration_error_gradient_difference(V,U);
        I=2-I; t=2-t;
    case 'gc'
        [I,t]=registration_error_gradient_correlation(V,U);
        I=1-I; t=1-t;
    case 'cc'
        [I,t]=registration_error_normalized_cross_correlation(V,U);
        I=1-I; t=1-t;
    case 'pi'
        [I,t]=registration_error_pattern_intensity(V,U);
        I=1-I; t=1-t;
    case 'ld'   
        [I,t]=registration_error_log_absolute_distance(V,U);
    otherwise
        error('Unknown error type')
end    
if(isnan(t)), warning('NaN in error image'); end

function [I,t]=registration_error_log_absolute_distance(V,U)
I=log(abs(V-U)+1);
t=sum(I(:))/numel(V);

function [I,t]=registration_error_normalized_cross_correlation(V,U)
Vvar=V-mean(V(:)); Uvar=U-mean(U(:));
I=(Vvar.*Uvar)/(sqrt(sum(Vvar(:).^2))*sqrt(sum(Uvar(:).^2)));
t=sum(I(:))/numel(V);

function [I,t]=registration_error_gradient_correlation(V,U)
if(ndims(U)==2)
    [Gx,Gy]=sobel2();
    I=(1/2)*(registration_error_normalized_cross_correlation(conv2(V,Gx,'same'),conv2(U,Gx,'same'))...
            +registration_error_normalized_cross_correlation(conv2(V,Gy,'same'),conv2(U,Gy,'same')));
else
    [Gx,Gy,Gz]=sobel3();
    I=(1/3)*(registration_error_normalized_cross_correlation(convn(V,Gx,'same'),convn(U,Gx,'same'))...
            +registration_error_normalized_cross_correlation(convn(V,Gy,'same'),convn(U,Gy,'same'))...
            +registration_error_normalized_cross_correlation(convn(V,Gz,'same'),convn(U,Gz,'same')));
end
t=sum(I(:))/numel(V);

function [Gx,Gy]=sobel2()
    Gx=[1 0 -1;2 0 -2;1 0 -1];  Gy=[1 2 1;0 0 0;-1 -2 -1];

function [Gx,Gy,Gz]=sobel3()
    Gx=zeros(3,3,3);Gy=zeros(3,3,3); Gz=zeros(3,3,3);
    Gx(:,:,1)=[-1 -3 -1;-3 -6 -3;-1 -3 -1]; Gx(:,:,2)=[ 0  0  0; 0  0  0; 0  0  0]; Gx(:,:,3)=[ 1  3  1; 3  6  3; 1  3  1];
    Gy(:,:,1)=[ 1  3  1; 0  0  0;-1 -3 -1]; Gy(:,:,2)=[ 3  6  3; 0  0  0;-3 -6 -3]; Gy(:,:,3)=[ 1  3  1; 0  0  0;-1 -3 -1];
    Gz(:,:,1)=[-1  0  1;-3  0  3;-1  0  1]; Gz(:,:,2)=[-3  0  3;-6  0  6;-3  0  3]; Gz(:,:,3)=[-1  0  1;-3  0  3;-1  0  1];    

function [I,t]=registration_error_squared_differences(V,U)
I=(V-U).^2;
t=sum(I(:))/numel(V);

function [I,t]=registration_error_differences(V,U)
I=(V-U);
t=sum(I(:))/numel(V);

function [I,t]=registration_error_pattern_intensity(V,U)
Idiff=V./mean(V(:))-U./mean(U(:));
o=0.3; %determines if grey-value varion is a structure (must be laster than noise)
r=5; numr=0; listr=[];
if(ndims(U)==2)
    for u=-r:r
        for v=-r:r
            if((u^2+v^2)<=r^2), numr=numr+1; listr(numr,:)=[u v]; end
        end
    end
    for u=-r:r
        for v=-r:r
            if((u^2+v^2)<=r^2), numr=numr+1; listr(numr,:)=[u v]; end
        end
    end
    SP1=zeros(size(U)-2*r,class(Idiff));
    for i=1:size(listr,1),
        u=listr(i,1); v=listr(i,2);
        SP1=SP1+o^2./(o^2+(Idiff(1+r:end-r,1+r:end-r)-Idiff(1+r+u:end-r+u,1+r+v:end-r+v)).^2);
    end
else
    for u=-r:r,
        for v=-r:r
            for w=-r:r
                if((u^2+v^2+w^2)<=r^2), numr=numr+1; listr(numr,:)=[u v w]; end
            end
        end
    end
    SP1=zeros(size(U)-2*r,class(Idiff));
    for i=1:size(listr,1),
        u=listr(i,1); v=listr(i,2); w=listr(i,3);
        SP1=SP1+o^2./(o^2+(Idiff(1+r:end-r,1+r:end-r,1+r:end-r)-Idiff(1+r+u:end-r+u,1+r+v:end-r+v,1+r+w:end-r+w)).^2);
    end
end
I=SP1./size(listr,1);
t=sum(I(:))/numel(V);


function [I,t]=registration_error_gradient_difference(V,U)
if(ndims(U)==2)
    a=mean(V(:))/mean(U(:));
    [Gx,Gy]=sobel2();
    Idiffv=conv2(V,Gx,'same')-a*conv2(U,Gx,'same');
    Idiffh=conv2(V,Gy,'same')-a*conv2(U,Gy,'same');
    Av=var(Idiffv(:));
    Ah=var(Idiffh(:));
    Sgdiff=(Av./(Av+Idiffv.^2))+(Ah./(Ah+Idiffh.^2));
else
    a=mean(V(:))/mean(U(:));
    [Gx,Gy,Gz]=sobel3();
    Idiffv=convn(V,Gx,'same')-a*convn(U,Gx,'same');
    Idiffh=convn(V,Gy,'same')-a*convn(U,Gy,'same');
    Idiffz=convn(V,Gz,'same')-a*convn(U,Gz,'same');
    Av=var(Idiffv(:));
    Ah=var(Idiffh(:));
    Az=var(Idiffz(:));
    Sgdiff=(Av./(Av+Idiffv.^2))+(Ah./(Ah+Idiffh.^2))+(Az./(Az+Idiffz.^2));    
end
I=Sgdiff;
t=sum(I(:))/numel(V);

function [I,t]=registration_error_mutual_info(V,U)
% This function t=registration_error_mutual_info(V,U) gives a registration error
% value based on mutual information (H(A) + H(B)) / H(A,B)

% Make a joint image histogram and single image histograms
bins=numel(V)^(1/ndims(V));
%bins=100; 

range=getrangefromclass(V);
if(isa(V,'double'))
    [hist12, hist1, hist2]=mutual_histogram_double(double(V),double(U),double(range(1)),double(range(2)),double(bins));
else
    [hist12, hist1, hist2]=mutual_histogram_single(single(V),single(U),single(range(1)),single(range(2)),single(bins));
end

% Calculate probabilities
p1=hist1./numel(V);
p2=hist2./numel(V);
p12=hist12./numel(V);

p1log=p1 .* log(p1+eps);
p2log=p2 .* log(p2+eps);
p12log=p12.* log(p12+eps);

% Calculate amount of Information
HA = -sum(p1log);
HB = -sum(p2log);
HAB = -sum(p12log(:));

% Studholme, Normalized mutual information
t=2-(HA+HB)/HAB;
I=[];

