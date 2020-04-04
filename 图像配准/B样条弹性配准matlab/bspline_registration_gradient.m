function [O_error, O_grad]=bspline_registration_gradient(O_grid,sizes,Spacing,I1,I2,options)
% Function registration_gradient. This function will calculate a registration
% error value and gradient after b-spline non-rigid registration
% of two images / volumes.
%
% [error, errorgrad]=bspline_registration_gradient(grid,sizes,Spacing,I1,I2,options)
%
% inputs,
%   grid: b-spline control grid created by make_init_grid.m, (and reshaped
%         to one long vector.)
%   sizes: sizes need tot reshape the grid to matrix format.
%   Spacing: The spacing in x,y (and z) direction of the b-spline grid
%           knots
%   I1 and I2: The input images, of which I1 is transformed.
%   options: Struct with options
%       ->type: Type of image similarity(error) measure used
%               (see image_difference.m) (default 'sd')
%       ->penaltypercentage Percentage of penalty smoothness (bending energy
%               of thin sheet of metal) (default 0.01)
%       ->step: Delta step used for error gradient (default 0.01)
%       ->centralgrad: Calculate derivatives with central instead of forward
%               gradient (default true)
% outputs,
%   error: The registration error value
%   errorgrad: The registration error gradient
%
% example,
%   I1=im2double(imread('lenag1.png')); 
%   I2=im2double(imread('lenag2.png'));
%   O_trans=make_init_grid([32 32],size(I1)); sizes=size(O_trans);
%   O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,[32 32],I1,I2,'sd'), O_trans(:), ...
%               optimset('Display','iter','GradObj','on','MaxIter',20,'DiffMinChange',0.1));
%   Icor=bspline_transform(reshape(O_trans,sizes),I1,Spacing); 
%   figure, imshow(I1), figure, imshow(I2), figure, imshow(Icor)
%
%
% This function is written by D.Kroon University of Twente (October 2008)

% Check/set input options
defaultoptions=struct('type','sd', 'penaltypercentage',0.01,'step',0.01,'centralgrad',true);
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('BsplineRegistrationGradient:unknownoption','unknown options found');
    end
end

% Type of image similarity(error) measure used
type=options.type;
% Percentage of penalty smoothness (bending energy of thin sheet of metal)
penaltypercentage=options.penaltypercentage;
% Delta step used for error gradient
step=options.step;
% Central gradient or (faster) forward gradient
centralgrad=options.centralgrad;


% Convert Grid vector to grid matrix
O_grid=reshape(O_grid,sizes);
            
% Transform the image with the B-spline grid
I_init=bspline_transform(O_grid,I1,Spacing);
   
% Calculate the current registration error
t=image_difference(I_init,I2,type);
O_error=t;

if(penaltypercentage>0)
    % Calculate penalty smoothness (bending energy of thin sheet of metal)
    if ( nargout > 1 )
        % If penalty gradient needed, also determine the gradient.
        [SO_error, SO_grad]=penalties_smoothness(O_grid,size(I1));
    else
        SO_error=penalties_smoothness(O_grid,size(I1));
    end

    % Add penalty to total error
    O_error=O_error+SO_error*penaltypercentage;
end

%
% Code below is only needed, when the error gradient is asked by the optimizer
%

% If gradient needed, also determine the gradient.
if ( nargout > 1 )
    if(ndims(I2)==2)
        O_grad=zeros(size(O_grid));
        for zi=0:3,
            for zj=0:3,
                % The variables which will contain the controlpoints for
                % determining a central registration error gradient
                O_gradpx=O_grid; O_gradpy=O_grid;
                if(centralgrad)
                    O_gradmx=O_grid; O_gradmy=O_grid;
                end

                %Set grid movements of every fourth grid node.
                for i=(1+zi):4:size(O_grid,1),
                    for j=(1+zj):4:size(O_grid,2),
                        O_gradpx(i,j,1)=O_gradpx(i,j,1)+step;
                        O_gradpy(i,j,2)=O_gradpy(i,j,2)+step;
                        if(centralgrad)
                            O_gradmx(i,j,1)=O_gradmx(i,j,1)-step;
                            O_gradmy(i,j,2)=O_gradmy(i,j,2)-step;
                        end
                    end
                end

                % Do the grid b-spline transformation for movement of nodes to
                % left right top and bottem.
                I_gradpx=bspline_transform(O_gradpx,I1,Spacing);
                I_gradpy=bspline_transform(O_gradpy,I1,Spacing);
                if(centralgrad)
                    I_gradmx=bspline_transform(O_gradmx,I1,Spacing);
                    I_gradmy=bspline_transform(O_gradmy,I1,Spacing);
                end
                
                for i=(1+zi):4:size(O_grid,1),
                    for j=(1+zj):4:size(O_grid,2),

                            % Calculate pixel region influenced by a grid node
                            irm=i-2; irp=i+2; jrm=j-2; jrp=j+2;
                            if(irm<1), irm=1; end
                            if(jrm<1), jrm=1; end
                            if(irp>size(O_grid,1)), irp=size(O_grid,1); end
                            if(jrp>size(O_grid,2)), jrp=size(O_grid,2); end

                            regAx=round(O_grid(irm,jrm,1)); regAy=round(O_grid(irm,jrm,2));
                            regBx=round(O_grid(irp,jrp,1)); regBy=round(O_grid(irp,jrp,2));

                            if(regAx<1), regAx=1; elseif(regAx>size(I1,1)), regAx=size(I1,1); end
                            if(regAy<1), regAy=1; elseif(regAy>size(I1,2)), regAy=size(I1,2); end
                            if(regBx<1), regBx=1; elseif(regBx>size(I1,1)), regBx=size(I1,1); end
                            if(regBy<1), regBy=1; elseif(regBy>size(I1,2)), regBy=size(I1,2); end
                            if(regAx>regBx), regAxt=regAx; regAx=regBx; regBx=regAxt; end
                            if(regAy>regBy), regAyt=regAy; regAy=regBy; regBy=regAyt; end
   
                            % Determine the registration error in the region
                            E_gradpx=image_difference(I_gradpx(regAx:regBx,regAy:regBy),I2(regAx:regBx,regAy:regBy),type);
                            E_gradpy=image_difference(I_gradpy(regAx:regBx,regAy:regBy),I2(regAx:regBx,regAy:regBy),type);
                            if(centralgrad)
                                E_gradmx=image_difference(I_gradmx(regAx:regBx,regAy:regBy),I2(regAx:regBx,regAy:regBy),type);
                                E_gradmy=image_difference(I_gradmy(regAx:regBx,regAy:regBy),I2(regAx:regBx,regAy:regBy),type);
                            else
                                E_grid=image_difference(I_init(regAx:regBx,regAy:regBy),I2(regAx:regBx,regAy:regBy),type);
                            end

                                % Calculate the error registration gradient.
                            if(centralgrad)
                                O_grad(i,j,1)=(E_gradpx-E_gradmx)/(step*2);
                                O_grad(i,j,2)=(E_gradpy-E_gradmy)/(step*2);
                            else
                                O_grad(i,j,1)=(E_gradpx-E_grid)/step;
                                O_grad(i,j,2)=(E_gradpy-E_grid)/step;
                            end
                    end
                end
            end
        end
        if(penaltypercentage>0)
            % Add smoothness penalty gradient
            O_grad=O_grad+SO_grad*penaltypercentage;
        end
        
        O_grad=O_grad(:);
    else
        O_grad=zeros(size(O_grid));
        for zi=0:3,
            for zj=0:3,
                for zk=0:3,
                    O_gradpx=O_grid; O_gradpy=O_grid; O_gradpz=O_grid;
                    if(centralgrad)
                        O_gradmx=O_grid; O_gradmy=O_grid; O_gradmz=O_grid;
                    end
    
                    %Set grid movements of every fourth grid node.
                    for i=(1+zi):4:size(O_grid,1),
                        for j=(1+zj):4:size(O_grid,2),
                            for k=(1+zk):4:size(O_grid,3),
                                O_gradpx(i,j,k,1)=O_gradpx(i,j,k,1)+step;
                                O_gradpy(i,j,k,2)=O_gradpy(i,j,k,2)+step;
                                O_gradpz(i,j,k,3)=O_gradpz(i,j,k,3)+step;
                                if(centralgrad)
                                    O_gradmx(i,j,k,1)=O_gradmx(i,j,k,1)-step;
                                    O_gradmy(i,j,k,2)=O_gradmy(i,j,k,2)-step;
                                    O_gradmz(i,j,k,3)=O_gradmz(i,j,k,3)-step;
                                end
                            end
                        end
                    end

                    % Do the grid b-spline transformation for movement of nodes to
                    % left right top and bottem.
                    I_gradpx=bspline_transform(O_gradpx,I1,Spacing);
                    I_gradpy=bspline_transform(O_gradpy,I1,Spacing);
                    I_gradpz=bspline_transform(O_gradpz,I1,Spacing);
                    if(centralgrad)
                        I_gradmx=bspline_transform(O_gradmx,I1,Spacing);
                        I_gradmy=bspline_transform(O_gradmy,I1,Spacing);
                        I_gradmz=bspline_transform(O_gradmz,I1,Spacing);
                    end

                    for i=(1+zi):4:size(O_grid,1),
                        for j=(1+zj):4:size(O_grid,2),
                            for k=(1+zk):4:size(O_grid,3),
                                % Calculate pixel region influenced by a grid node
                                irm=i-2; irp=i+2; 
                                jrm=j-2; jrp=j+2; 
                                krm=k-2; krp=k+2;
                                if(irm<1), irm=1; end
                                if(jrm<1), jrm=1; end
                                if(krm<1), krm=1; end
                                if(irp>size(O_grid,1)), irp=size(O_grid,1); end
                                if(jrp>size(O_grid,2)), jrp=size(O_grid,2); end
                                if(krp>size(O_grid,3)), krp=size(O_grid,3); end

                                regAx=round(O_grid(irm,jrm,krm,1)); regAy=round(O_grid(irm,jrm,krm,2)); regAz=round(O_grid(irm,jrm,krm,3));
                                regBx=round(O_grid(irp,jrp,krp,1)); regBy=round(O_grid(irp,jrp,krp,2)); regBz=round(O_grid(irp,jrp,krp,3));

                                if(regAx<1), regAx=1; elseif(regAx>size(I1,1)), regAx=size(I1,1); end
                                if(regAy<1), regAy=1; elseif(regAy>size(I1,2)), regAy=size(I1,2); end
                                if(regAz<1), regAz=1; elseif(regAz>size(I1,3)), regAz=size(I1,3); end
                                if(regBx<1), regBx=1; elseif(regBx>size(I1,1)), regBx=size(I1,1); end
                                if(regBy<1), regBy=1; elseif(regBy>size(I1,2)), regBy=size(I1,2); end
                                if(regBz<1), regBz=1; elseif(regBz>size(I1,3)), regBz=size(I1,3); end
                                if(regAx>regBx), regAxt=regAx; regAx=regBx; regBx=regAxt; end
                                if(regAy>regBy), regAyt=regAy; regAy=regBy; regBy=regAyt; end
                                if(regAz>regBz), regAzt=regAz; regAz=regBz; regBz=regAzt; end
                            
                                % Determine the registration error in the region
                                E_gradpx=image_difference(I_gradpx(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                E_gradpy=image_difference(I_gradpy(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                E_gradpz=image_difference(I_gradpz(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                if(centralgrad)
                                    E_gradmx=image_difference(I_gradmx(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                    E_gradmy=image_difference(I_gradmy(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                    E_gradmz=image_difference(I_gradmz(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                else
                                    E_grid=image_difference(I_init(regAx:regBx,regAy:regBy,regAz:regBz),I2(regAx:regBx,regAy:regBy,regAz:regBz),type);
                                end
                                
                                % Calculate the error registration gradient.
                                if(centralgrad)
                                    O_grad(i,j,k,1)=(E_gradpx-E_gradmx)/(step*2);
                                    O_grad(i,j,k,2)=(E_gradpy-E_gradmy)/(step*2);
                                    O_grad(i,j,k,3)=(E_gradpz-E_gradmz)/(step*2);
                                else
                                    O_grad(i,j,k,1)=(E_gradpx-E_grid)/step;
                                    O_grad(i,j,k,2)=(E_gradpy-E_grid)/step;
                                    O_grad(i,j,k,3)=(E_gradpz-E_grid)/step;
                                end
                            end
                        end
                    end
                end
            end
        end
        if(penaltypercentage>0)
            % Add smoothness penalty gradient
            O_grad=O_grad+SO_grad*penaltypercentage;
        end
        
        O_grad=O_grad(:);
    end
end   



    
    