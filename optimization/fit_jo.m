function [output] = fit_jo(funName, params, freeList, varargin)
% [params,err] = fit(funName, params, freeList, var1, var2,...)

%% this is a modified version of the 'fit' function, which is described below.
% modified by J. Carpenter, 2020. > @todo: describe changes made
%
% Helpful interface to MATLAB's 'fminsearch' function.
%
% Inputs:
%   funName        Function to be optimized. Must have form
%                  [err] = <funName>(params, var1, var2, ...)
%
%   params         A structure of with field names that correspond to 
%                  parameter values for fitted function. Params are allowed
%                  to be matrices.
%       options    A structure with options for MATLAB's fminsearch program
%                  (see OPTIMSET)
%
%   freeList       Cell array containing list of parameter names (strings)
%                  to be free in fitting. Free strings can contain certain
%                  values / ranges within the 'params' matrices. For 
%                  example, the following are valid.
%
%                  {'x(1)','y(3:4)', 'z(1:2,4:5)'}
%
%   var<n>         Extra variables to be sent into fitted function
%                  'funName'
%
% Outputs:
%   params         A structure with best fitting parameters as fieldnames
%
%   err            Error value at minimum, numeric
%
% Notes:
% - Dependencies: params2var.m, var2params.m, fitFunction.m

% Written by Geoffrey M. Boynton, Summer of '00
% Edited by Kelly Chang, February 10, 2017

%% Input Control

options = []; % options for fminsearch (see OPTIMSET)
if isfield(params, 'options')
    options = params.options;
end

if isempty(freeList)
    freeList = fieldnames(params);
end

%% Fit Function and Calculate Final Error

% turn free parameters in to 'var'
vars = params2var(params, freeList);

% calling fminsearch
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
options = optimset('Display','off');
[vars,fval,exitflag,out, save] = fminsearch_mod('fitFunction', vars, options, funName, params, freeList, varargin);
% [results,vars,fval,exitflag,out] = evalc("fminsearch('fitFunction', vars, options, funName, params, freeList, varargin)");

% if you want to get the data from the plot
% figHandles = findobj('Type', 'figure');
% ax = gca;
% fig = get(ax,'children');

% assign final parameters into 'params'
params = var2params(vars, params, freeList);

% evaluate the function 'funName' for error at minimum
err = fitFunction(vars, funName, params, freeList, varargin);

% make an output struct to store everything in
% output.plot.X = fig.XData;
% output.plot.Y = fig.YData;
% output.plot.title = ax.Title.String;
output.params = params;
output.error = err;
output.exitflag = exitflag;
output.fval = fval;
output.out = out;
output.saved = save;
% output.out.results = results;
end

%% SCRATCH
% [results, vars,fval,exitflag,output] = evalc("fminsearch('fitFunction', vars, options, funName, params, freeList, varargin)");
% [results, vars,fval,exitflag,output] = evalc("fminsearch_mod('fitFunction', vars, options, funName, params, freeList, varargin)");
% options = optimset('Display','final','PlotFcns',@optimplotfval);
% [vars, ~, history] = fminsearch_outfun('fitFunction', vars, options, funName, params, freeList, X, Y, H, R);

% history = [];
% options = optimset('OutputFcn', @myoutput);
% [vars,fval,exitflag,output] = fminsearch('fitFunction', vars, options, funName, params, freeList, varargin);
%     

% options = optimset('OutputFcn', @myoutput);
% [vars,fval,exitflag,output] = fminsearch_mod('fitFunction', vars, options, funName, params, freeList, varargin);

