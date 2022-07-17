function [r t_max] = ref(t)
% Defines reference signal to be tracked by muscles w/ Pontryagin control,
% relative to the minimum force produced by the muscles due to operation
% within the linear region
% Default reference is overwritten if global anonymous function refFun is
% defined
% r [N] reference force at time t, relative to minimum force produced
% t_max [s] duration of tracked reference

t_max = 200;

global setRef
if isempty(setRef); setRef=1; end
% setRef==1 is default reference
if setRef == 2
    %t_max = 80;
    f_ref = 0.05; % Hz square-wave with offset
    r = 1*(mod(t+1/(4*f_ref),1/f_ref)>1/(2*f_ref));
elseif setRef == 3
    %t_max = 100;
    N = 5;
    f = linspace(0.01, 0.1, N);
    r = sum(sin(2*pi*f.'*(t+81) + 0.9058...
       - (pi*[1:N].^2/N).'*ones(1,length(t))))/5.8579 + 0.5024;
elseif setRef == 4
    f = 0.05; % Hz
    t_max = 2/f; %overwrite t_max
    r = 0.5*(-cos(2*pi*f*t)/2+0.5);
elseif setRef == 5
    f = 0.1; % Hz
    t_max = 2/f; %overwrite t_max
    r = 0.5*(-cos(2*pi*f*t)/2+0.5);
else %setRef == 1
    f_ref = 0.02; % Hz square-wave with offset
    r = 1*(mod(t+1/(4*f_ref),1/f_ref)>1/(2*f_ref));
end

end