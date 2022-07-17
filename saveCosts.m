clc
clear all
close all

% See "Redundant Actuation of Twisted and Coiled Polymer Muscles to Improve
% Tracking Performance" MSc Thesis by Ian McGinty for context. This thesis
% is available through the TU Delft Library.

% The following matlab code iteratively halves weighting matrix scalar R0
% until tracking cost converges between 2 successive iterations within 2%
% Scalar R0 dictates the relative contributions of tracking cost and
% control effort cost to total cost, which optimal control minimizes.
% Using a provided guess mesh and solution, getOptimalControl solves the
% split-boundary-value problem defined the thesis, which os a function of
% R0. The code is robust to mesh density limit errors. If mesh density
% limit is exceeded, the initial step in R0 is broken into recursively
% smaller subintervals.

% Define TCPA configruation charactersisitcs:
ant = false; % bundle of n unidirectional or antagonsitic pairs of TCPMs
ref = 1; % reference signal defined by ref.m, enumerated per Fig 9
d = [0.0004 0.0004]; % [m] nx1 vector of fiber diameters
C = [2 2]; % nx1 vector of spring indices
n = length(d);
% See initializeTCPAs.m to define TCPMs with different twist insertions,
% lengths, etc.




global setRef; setRef = ref;

if ant;     file = 'ant_';
else        file = 'uni_'; end
file = strcat(file,num2str(setRef),'_n');
file = strcat('./costs/',file,num2str(n),'costs.mat');
disp(strcat('saved to:',file));

% Convergance critertia with respect to R
j_max = 7; % max number of times R is lowered
costTol = 0.02; % convergance between values of R
R0 = 0.2;
if ref==3
    R0 = 0.1;
    j_max = 25;
end
if n == 4
    R0 = 0.4;
end
relTol = costTol/10; % convergance within each value of R

output.costTol = costTol;
output.relTol = relTol;
output.R0 = R0./(2.^[0:j_max-1]);
output.warningFlags = 0;

output.d = d;
output.C = C;

x = [d C];

tic

i = 1;
initializeTCPAs(x,n);
setR = R0;

% Initial solution for R = R0
% Use default guess for bvp4c
[~, ~, ~, costE, costT, ~, ~, sol] = ...
    getOptimalControl(R0,relTol,ant,false,[]);
costError(i,1) = costE;
costTotal(i,1) = costT;

disp(strcat('R=',num2str(setR),' for x=',num2str(x)))

% Iteratively lower R until perf cost converges
j = 2;
lowerR = true;
tryHalfStep = true;
while lowerR % continue until convergance or error
    setR = R0/(2^(j-1));
    disp(strcat('R=',num2str(setR),' for x=',num2str(x)))

    % lowers R by a factor of 2
    % if mesh density limit hit, iteratively lowers R by a factor of
    % 1.5 (or 1.25 or 1.125, etc) to ensure mesh density limit is not
    % hit
    [costE, costT, sol, ~] = recursivelyLowerR(0, 1, setR, relTol, ant, sol);
    costError(i,j) = costE;
    costTotal(i,j) = costT;

    relCostChange = (costE - costError(i,j-1))/costE;
    j = j+1;

    warnMsg = lastwarn;
    if ~isempty(warnMsg) % bvp4c still hits mesh density limit
        lowerR = false;
        costError(i,j-1:j_max) = NaN;
        costTotal(i,j-1:j_max) = NaN;
        msg{i,1} = strcat('Convergence failed: bvp4c warning. Last entry removed. Converged to:',...
            num2str(abs(costError(i,j-3)/costError(i,j-2)-1)*100),'%');
        warning(msg{i,1})
        output.warningFlags = output.warningFlags + 1;
    elseif (abs(relCostChange) < costTol) % sol converged w.r.t. R
        lowerR = false;
        costError(i,j:j_max) = NaN;
        costTotal(i,j:j_max) = NaN;
        msg{i,1} = strcat('Converged to:',...
            num2str(abs(relCostChange)*100),'%');
    elseif j > j_max % max number of iterations exceeded
        lowerR = false;
        msg{i,1} = strcat('Convergence failed: R < R_min. Converged to:',...
            num2str(abs(relCostChange)*100),'%');
        warning(msg{i,1});
        output.warningFlags = output.warningFlags + 1;
    end

    output.costError = costError;
    output.costTotal = costTotal;
    save(file,'output');
end
% Save after every geometry
output.msg = msg;
save(file,'output');

% Prints graph (or outputs) results for the final value of R0
getOptimalControl(setR,relTol,ant,true,sol);

toc


%%
function [costE, costT, sol, x] = recursivelyLowerR(x, dx, setR, relTol, ant, sol)
% increase x until x=1 such that 2*setR/(1+x) = setR
% terminates if recursive function sets x to NaN
while x < 1 && ~isnan(x)
    if dx ~= 1; disp(strcat('R=',num2str(2*setR/(1+x+dx)))); end
	if dx >= 1/128 % continue if step is not too small
        warning('')
        [~, ~, ~, costE, costT, ~, ~, newSol] = ...
                getOptimalControl(2*setR/(1+x+dx),relTol,ant,false,sol);
        warnMsg = lastwarn;
        if isempty(warnMsg) % if convergence successful, step x
            x = x + dx;
            sol = newSol;
        else % if convergence unsuccessful, decrease step size
            [costE, costT, sol, x] = recursivelyLowerR(x, dx/2, setR, relTol, ant, sol);
            x = x + dx/2;
        end
    else % terminates if step is too small
        costE = NaN;
        costT = NaN;
        x = NaN;
    end
end
end
