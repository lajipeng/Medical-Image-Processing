function [x,fval,exitflag,output,grad] =fminsd(func,x_init,optim)
%FMINSD finds a local minimum of a function of several variables,
%   using Steepest Gradient Descent optimization.
%   X = FMINSD(FUN,X0) starts at X0 and attempts to find a local minimizer
%   X of the function FUN. FUN accepts input X and returns a scalar
%   function value F evaluated at X. X0 can be a scalar, vector or matrix. 
%
%   The speed of this optimizer can be improved by also providing
%   the gradient at X. Write the FUN function as follows
%   function [f,g]=FUN(X)
%       f , value calculation at X;
%   if ( nargout > 1 )
%       g , gradient calculation at X;
%   end
%
%   X = FMINSD(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details, you 
%   can also use a normal STRUCT when no optimalization toolbox is
%   available. 
%   Used options are Display, TolX, TolFun, GradObj, MaxIter, MaxFunEvals, 
%   DiffMaxChange, DiffMinChange and OutputFcn.
%
%   [X,FVAL] = FMINSD(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FMINSD(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of FMINUNC. Possible values of EXITFLAG 
%   and the corresponding exit conditions are
%
%     1  Magnitude of gradient smaller than the specified tolerance. 
%     2  Change in X smaller than the specified tolerance.
%     3  Change in the objective function value smaller than the specified 
%         tolerance (only occurs in the large-scale method).
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Algorithm terminated by the output function.
%    -2  Line search cannot find an acceptable point along the current
%         search direction (only occurs in the medium-scale method).
%   
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSD(FUN,X0,...) returns a structure 
%   OUTPUT with the number of iterations taken in output.iteration, the 
%   number of function evaluations in OUTPUT.funcCount, the algorithm used 
%   in OUTPUT.algorithm, 
%
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINSD(FUN,X0,...) returns the value 
%   of the gradient of FUN at the solution X.
%
%
%   Examples
%     FUN can be specified using @:
%        X = fminsd(@myfun,2)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x) + 3;
%
%     To minimize this function with the gradient provided, modify
%     the function myfun so the gradient is the second output argument:
%        function [f,g] = myfun(x)
%         f = sin(x) + 3;
%         g = cos(x);
%     and indicate the gradient value is available by creating an options
%     structure with OPTIONS.GradObj set to 'on' (using OPTIMSET):
%        options = optimset('GradObj','on');
%        x = fminsd(@myfun,4,options);
%
%   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, FMINUNC, @, INLINE.
%
%   Function is written by D.Kroon University of Twente (October 2008)

% Read Optimalisation Parameters
defaultopt = struct('Display','final','TolX',1e-6,'TolFun',1e-6,'GradObj','off','MaxIter',400,'MaxFunEvals',100*numel(x_init),'DiffMaxChange',1e-1,'DiffMinChange',1e-8,'OutputFcn',[]);
if (~exist('optim','var')) 
    optim=defaultopt;
else
    f = fieldnames(defaultopt);
    for i=1:length(f),
        if (~isfield(optim,f{i})||(isempty(optim.(f{i})))), optim.(f{i})=defaultopt.(f{i}); end
    end
end

% Make a long vector from x_init but store the dimensions for 
% function calls
xsizes=size(x_init);
x=x_init(:);

% Initialize variables
output.iteration=0;
output.funcCount=0;
output.stepsize=1;
output.message='';
output.algorithm='Steepest Gradient Descent';
output.fval=0;
output.gradient=0;

% Calculate starting error and gradient
output=gradient_function(func,x,output,optim,xsizes);

% Execute Output function
if(~isempty(optim.OutputFcn))
    stopt=optim.OutputFcn(x, output, 'init');
    if(stopt), exitflag=-1; end
end

if(strcmp(optim.Display,'iter'))
    disp('    Iteration  Func-count      f(x)      Step-size');
    disp(['        '  int2str(output.iteration) '          ' int2str(output.funcCount) '          ' num2str(output.fval) '       ' ]);
end

output=init_stepsize(func,x,output,xsizes);

grad_old=output.gradient; x_old=x; fval_old=output.fval; 
while(true)
    % Find a good step size in the direction of the gradient (linesearch)
    [output,exitflag]=linesearch(func,x_old,grad_old,fval_old,output,optim,xsizes);
  
    % Update x with step in direction of gradient
    x=x_old-output.stepsize*grad_old(:)/sqrt(sum(grad_old(:).^2));
    
    % Calculate the new error and new gradient
    [output]=gradient_function(func,x,output,optim,xsizes);

    % Set exit flags if difference in error or x between iterations is too small
    if(sum(abs(output.gradient(:)))<optim.TolFun), exitflag=1; end
    if(abs(sum(x.^2)-sum(x_old.^2))<optim.TolX), exitflag=2; end
            
    % Show the current iteration
    output.iteration=output.iteration+1; 
    if(strcmp(optim.Display,'iter'))
        disp(['        '  int2str(output.iteration) '          ' int2str(output.funcCount)  '          '  num2str(output.fval)  '       '  num2str(output.stepsize)]);
    end
    
    % Execute Output function
    if(~isempty(optim.OutputFcn))
        stopt=optim.OutputFcn(x, output, 'iter');
        if(stopt), exitflag=-1; end
    end
    
    % Set exit flags if number of iterations or Function eval is too large
    if(output.iteration>=optim.MaxIter), exitflag=0; end
    if(output.funcCount>=optim.MaxFunEvals), exitflag=0; end

    % Check if exitflag is set
    if(~isempty(exitflag)), break, end;
    
    % Keep the variables for next iteration
    x_old=x; grad_old=output.gradient; fval_old=output.fval;
end

fval=output.fval;
grad=output.gradient;

% Reshape x to original shape
x=reshape(x,xsizes);

% get exit message
output.message=getexitmessage(exitflag);

% Execute Output function
if(~isempty(optim.OutputFcn))
    optim.OutputFcn(x, output, 'done');
end
    
% Display final results
if(~strcmp(optim.Display,'off'))
    disp('    Optimizer Results')
    disp(['        Exit message : ' output.message]);
    disp(['        iterations : '  int2str(output.iteration)]);
    disp(['        Function Count : ' int2str(output.funcCount)]);
    disp(['        Minimum found : ' num2str(fval)]);
    disp(['        Minimum stepsize : ' num2str(output.stepsize)]);
end

function message=getexitmessage(exitflag)
    switch(exitflag)
        case 1, message='Magnitude of gradient smaller than the specified tolerance TolFun.';
        case 2, message='Change in x was smaller than the specified tolerance TolX.'; 
        case 3, message='Change in the objective function value was less than the specified tolerance';
        case 0, message='Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
        case -1, message='Algorithm was terminated by the output function.';
        case -2, message='Line search cannot find an acceptable point along the current search';
        otherwise, message='Undefined exit code';
    end

function [output,exitflag]=linesearch(func,x_old,grad_old,fval_old,output,optim,xsizes)
    exitflag=[]; itw=0; found=false; stepsizeb=0.01; stepsize=output.stepsize;
    while((~found)||itw<7)
        % Move the grid points with use of gridpoint gradient
        x=x_old-stepsize*grad_old(:)/sqrt(sum(grad_old(:).^2));

        % Calculate the error registration gradient of the gridpoints
        [fval]=func(reshape(x,xsizes));  output.funcCount=output.funcCount+1; if(output.funcCount>=optim.MaxFunEvals), exitflag=0; break; end
                
        % Update step value
        if(fval_old<fval), stepsize=stepsize*0.6;  end
        if(fval_old>fval), stepsizeb=stepsize; found=true; fval_old=fval; stepsize=stepsize*1.15;   end
                
        % Update number of loop iterations
        itw=itw+1; if(itw>100), exitflag=-2; break; end 
   end
   output.stepsize=stepsizeb;
    
            
function output=gradient_function(func,x,output,optim,xsizes)
    % Call the error function for error and gradient
    if(strcmp(optim.GradObj,'on'))
        [fval, grad]=func(reshape(x,xsizes)); output.funcCount=output.funcCount+1;
    else
        % Calculate gradient if not provided by the function
        grad=zeros(length(x),1);
        [fval]=func(reshape(x,xsizes));
        gstep=output.stepsize/1000; 
        if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
        if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
        for i=1:length(x),
            x_temp=x; x_temp(i)=x_temp(i)+gstep;
            [fval_g]=func(reshape(x_temp,xsizes));  output.funcCount=output.funcCount+1;
            grad(i)=(fval_g-fval)/gstep;
        end
    end
    %grad=reshape(grad,[12 12 2]); grad(1:2,:,:)=0; grad(:,1:2,:)=0; grad(11:12,:,:)=0; grad(:,11:12,:)=0;
    output.fval=fval;
    output.gradient=reshape(grad,xsizes);
    

function output=init_stepsize(func,x_init,output,xsizes)
    % Find best stepsize in the range of 0.0001 to 1000
    fval_old=inf; output.stepsize=0; grad_init=output.gradient;
    for i=-4:3
        stept=10^i;
        x=x_init-stept*grad_init(:)/sqrt(sum(grad_init(:).^2)); 
        [fval]=func(reshape(x,xsizes)); output.funcCount=output.funcCount+1;
        if(fval<fval_old), output.stepsize=stept; fval_old=fval; end
    end     
                                        