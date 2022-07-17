function [C_twist, a_z, a_t, a_r, a_zt] = getMaterialProperties(beta_f, T)
% returns twisted stiffness matrix and coefficients of thermal expansion
% based on the bias angle of the current layer
% beta_f fiber bias angle in radians, CCW is positive
% T is temperature is in ºC

persistent print;
if isempty(print)
    print = true;
end

if(print)
disp('>')
disp('> getMaterialProperties.m')
disp('> Temperature dependence of material properties not accounted for')
end
% 1 is the radial direction
% Transverse isotropy on the 2-3 plane
% 1-2-3 in the twisted bases corresponds to z-theta-r in fiber basis

if(print)
disp('> Only a1 is temp dependent, accurate up to 180ºC')
disp('> Using material properties from 66_Zhao')
end
E1 = 2.66E9; % [Pa] Youngs modulus is axial direction
E2 = 0.56E9; % [Pa] Youngs modulus is radial direction, ~20-25% E1 [66]
% v_ij corresponds to a contraction in j when an extension is applied in i
% v_ij/E_i = v_ji/E_j
% https://en.wikipedia.org/wiki/Poisson%27s_ratio
v12 = 0.48; % = v13, tensile Poisson ratio
v23 = 0.26; % = v32, transverse Poisson ratio
G12 = 1.13E9; % = G13 [Pa] axial shearing modulus
a1 = getAlpha1(T);
a2 = 8.1E-5;

%{
% "Modeling of TCP Muscles for Understanding Actuation Behavior" Tadesse
if(print)
disp('> Material properties with augmented with Tadesse temp dependence')
end
T_amb  = getT_amb();
To = 20;
E1 = E1*(T/(T_amb+To))^-1.079;
E2 = E2*(T/(T_amb+To))^-1.079;
G12 = G12*(T/(T_amb+To))^-1.079;
%}
% Compliance matrix, S: strain = S*stress 
% See wiki for compliant matrix of a transversly isotropic material:
%   https://en.wikipedia.org/wiki/Poisson%27s_ratio
% 5 independent elasitc matierials required: E1, E2, G12, v12, v23
G23 = E2/(2+2*v23);
S = diag([E1 E2 E2 G23 G12 G12].^-1);
S(1,2) = -v12/E1; S(2,1) = S(1,2);
S(1,3) = -v12/E1; S(3,1) = S(1,3);
S(2,3) = -v23/E2; S(3,2) = S(2,3);
C = inv(S);
%{
% "Mechanical Relaxations and Moduli of Oriented Nylon66 and Nylon 6"
% Choy 1984, Dry Nylon 66 Coeff at temps of -40, 40, 90, 150ºC
S11i = [0.120	0.182	0.278	0.82]*10^-9;
S12i = [-0.046	-0.07	-0.115	-0.36]*10^-9;
S22i = [0.212	0.290	0.41	0.97]*10^-9;
S23i = [-0.106	-0.151	-0.216	-0.54]*10^-9;
S55i = [0.67	0.9	1.15	3.02]*10^-9;

Ti = [-40 40 90 150];
S11 = 10.^interp1(Ti,log10(S11i),T);
S12 = 10.^interp1(Ti,log10(S12i),T);
S22 = 10.^interp1(Ti,log10(S22i),T);
S23 = 10.^interp1(Ti,log10(S23i),T);
S55 = 10.^interp1(Ti,log10(S55i),T);
S44 = 2*(S22-S23);
S = diag([S11 S22 S22 S44 S55 S55]);
S(1,2) = S12; S(2,1) = S(1,2);
S(1,3) = S12; S(3,1) = S(1,3);
S(2,3) = S23; S(3,2) = S(2,3);


T1 = [-44,-1,37,50,64,79,97,116,139,159];
alpha1 = [2.89,2.59,2.67,3.23,3.79,3.27,1.81,-0.86,-4.73,-8.65]*10^-5;
%a1 = interp1(T1,alpha1,T);
T2 = [-41,1,32,64,92,116,137,160];
alpha2 = [10.84,12.81,15.50,19.56,24.90,30.48,36.35,43.68]*10^-5;
%a2 = interp1(T2,alpha2,T);

if T>max(Ti) || T<min(Ti)
    warning('Temp out of range of elastic material properties defined by Choy 1984')
end

%C = inv(S);
%}
% Only looking at non-zero stress states z, theta, r, z-theta
% Removing rows and columns 5 and 6 does not affect matrix inversion
% or twist

% Return rotated stiffness matrix and coefficients of thermal expansion
s = sin(beta_f);
c = cos(beta_f);
% Counterclockwise rotation about 3 axis
R3 =[c^2    s^2 0   0   0   2*c*s;
     s^2    c^2 0   0   0   -2*c*s;
     0      0   1   0   0   0;
     0      0   0   c   s   0;
     0      0   0   -s  c   0;
     -c*s   c*s 0   0   0   c^2 - s^2];
C_twist = R3*C*R3.';
a_z = a1*c^2 + a2*s^2;
a_t = a1*s^2 + a2*c^2;
a_r = a2;
a_zt = 2*(a2-a1)*c*s; % Opposite sign as 66_Zhao b/c this is CCW rotation
% factor of 2 due to definition of engineering vs true strain

% Material property rotation overview
% http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
%{
clc
clear all

syms c s
syms C [3 3];
syms C4_4 C5_5 C6_6;
C(4:6, 4:6) = diag([C4_4 C5_5 C6_6]);
C(2,1) = C(1,2);
C(3,1) = C(1,3);
C(3,2) = C(2,3);

s = -s; % Switch from CCW to CW rotation

syms a11 a22 a12
a =  diag([a11 a22 a22]);
 
% Counterclockwise rotation about 3 axis
R3 =[c  -s  0
     s  c   0
     0  0   1];
a_twist = R3*a*R3.'; % This may supposed to be R3.'*a*R3...

% Counterclockwise rotation about 2 axis
R2 =[c^2    0   s^2 0   2*c*s   0;
     0      1   0   0   0       0;
     s^2    0   c^2 0   -2*c*s  0;
     0      0   0   c   0       -s;
     -c*s   0   c*s 0   c^2-s^2 0;
     0      0   0   s   0       c];
 
% Counterclockwise rotation about 3 axis
R3 =[c^2    s^2 0   0   0   2*c*s;
     s^2    c^2 0   0   0   -2*c*s;
     0      0   1   0   0   0;
     0      0   0   c   s   0;
     0      0   0   -s  c   0;
     -c*s   c*s 0   0   0   c^2 - s^2];
 
C_twist = R3*C*R3.';

% Following match 25_Li and 66_Zhao when s = -s (clockwise rotation)
% Note each 25 and 66 have 1 different incorrect term
expand(C_twist(1,1))
expand(C_twist(1,2))
expand(C_twist(1,3))
expand(C_twist(1,6))
expand(C_twist(2,2))
expand(C_twist(2,3))
expand(C_twist(2,6))
expand(C_twist(3,3))
expand(C_twist(3,6))
expand(C_twist(6,6))
% a11 and a22 match, but a22 off factor of 2 
% This indicates wrong rotation direction, and true strain factor of 2
% included in a_twist matrix, but not in compliance matrix; [25] and [66]
% are sloppy
a_twist

%}

print = false;

end