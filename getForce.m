function [F, dF, dL_th] = getForce(T,d,l0,l,L0,L,N,blockedTorque)
% Returns force F calculated using Castigliano's method
% T [ºC] temperature for if material props are temp dependent
% d [m] fiber diameter
% l0 [m] initial length of twisted fiber
% L0 [m] initial length of TCPA 
% L [m] stretched length of TCPA 
% N number of TCPA coils
% strain_z axial strain of the fiber such that l = l0*(1+strain_z)
% blockedTorque [Nm]

persistent print;
if isempty(print)
    print = true;
end

if(print)
disp('>')
disp('> getForce.m')
disp('> Averege twisted material properties are emperically determined!')
disp('> Change in fiber length l is accounted for; change in fiber OD d is not')
disp('> Fiber untwist of beta_f_d due to change in l and d is assumed negligible')
disp('> Stiffness contribution of heating wire is ignored')
end

% TCPA material properties
% Twisted material properties average over fiber cross section
% These may be temperature dependent
Ez_bar = 2.2E9; % [Pa] emperically determined in 66_Zhao
Gzt_bar = 0.43E9; % [Pa] emperically determined in 66_Zhao

sin_beta_c0 = L0/l0;
sin_beta_c = L/l;

% Castiglianos Method: dL = f11*F - f12*blockedTorque;
% Coil Geometry: dL = l*sin_beta_c - l0*sin_beta_c0

A = 8*l^3/(pi^3*d^4*N^2);
B = 4*l/(pi*d^2);
C = 16*l^2/(pi^2*d^4*N);

f11 = A*((1-sin_beta_c^2)^2/Gzt_bar + ...
    2*sin_beta_c^2*(1-sin_beta_c^2)/Ez_bar) + ...
    B*((1-sin_beta_c^2)/Gzt_bar + sin_beta_c^2/Ez_bar);
f12 = C*(1-sin_beta_c^2)/Gzt_bar;

dL = l*sin_beta_c - l0*sin_beta_c0;
dL_th = f12*blockedTorque; % contraction due to thermal force
dF = dL_th/f11;
F = dL/f11 + dF;

print = false;

end
