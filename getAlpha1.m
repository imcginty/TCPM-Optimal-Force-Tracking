function [a1] = getAlpha1(T)
% Returns coefficient of thermal expansion per 66_Zhao Fig 9
% T is temperature is in ºC

% Cos approximation valid until ~140 C
% a1 = (-0.2 + 3.05/2*(cos((T-50)*pi/95)-1))*10^-4;

persistent print;
if isempty(print)
    print = true;
end


if T > 180
error('getAlpha1.m: T > 180; alpha_1 not defined')
end
if T < 40 && T >= 20 && print
disp('>')
disp('> getAlpha1.m')
disp('> WARNING: T < 40; alpha_1 extrapolated below 40ºC')
print = false;
end
if T < 20
error('getAlpha1.m: T < 20; alpha_1 not defined')
end

temp = [20:10:180];
alpha = [   -0.25;
            -0.25;
            -0.25;
            -0.20;
            -0.28;
            -0.47;
            -0.81;
            -1.27;
            -1.77;
            -2.29;
            -2.75;
            -3.10;
            -3.25;
            -3.15;
            -2.8;
            -2.2;
            -1.25].'*10^-4;
a1 = interp1(temp,alpha,T);

end

