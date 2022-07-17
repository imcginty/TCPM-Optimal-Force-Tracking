function [blockedTorque, T, stress_zt, r, delta_z] = ...
    getBlockedTorque(beta_f_d,d,T0,T_max)
% Returns a vector of blocked torque [Nm] and associated temp T [ºC]
% Optionally returns shear stress [Pa] @ T=T_max as funtion of radius r
% Optionally returns hypothetical axial strain delta_z as function of T
% Note delta_z is independent of the sign of beta_f_d
% Material properties defined in getMaterialProperties
% beta_f_d [rad] fiber bias angle from vertical @ r=d/2, CCW is +
% beta_f_d sign defintion determines (+) homo vs (-) heterochiral 
% beta_f_d sign does not affect magnitude of stress_zt or blockedTorque
% d [m] fiber diameter
% T0, T_max [ºC] initial and maximum temperature of interest

persistent print;
if isempty(print)
    print = true;
end

if(print)
disp('>')
disp('> getBlockedTorque.m')
disp('> Method of 66_Zhao assumes average axial stress across fiber = 0')
end

% Integrate moment induced by shear stress from T0 to T_max
% Define integration steps in T
T = [logspace(log10(T0), log10(100),100) linspace(101, T_max,10)];
blockedTorque = 0;
delta_z = 0;

for j = 2:length(T)

% For each temp step, find moment induced by shear stress at temp T by
% integrating radially dependent moment induced by shear stress across
% fiber cross section
% Define integration steps in r
r = linspace(0,d/2,10);

% Per 66_Zhao, first calculate axial fiber displacement (Eq 11) at temp T
% getMatProp called twice; could be more efficient, but less elegant
for i = 1:length(r)
    % Find radially dependent material properties
    beta_f = atan(2*r(i)/d*tan(beta_f_d));
    [C, a_z, a_t, a_r, a_zt] = getMaterialProperties(beta_f, T(j));
    % Find integrand of Eq 11
    num(i) = r(i)*(C(1,1)*a_z + C(1,2)*a_t + C(1,6)*a_zt);
    den(i) = r(i)*C(1,1);
end
delta_z_dT(j) = trapz(r,num)/trapz(r,den);
delta_z(j) = ...
    trapz([T(j-1) T(j)], [delta_z_dT(j-1) delta_z_dT(j)]) + delta_z(j-1);

% Per 66_Zhao, then find radially dependent shear stress
for i = 1:length(r)
    %Find radially dependent material properties
    beta_f = atan(2*r(i)/d*tan(beta_f_d));
    [C, a_z, a_t, a_r, a_zt] = getMaterialProperties(beta_f, T(j));
    % See Eq 7
    stress_zt_dT(i,j) = C(1,6)*(delta_z_dT(j)-a_z) - C(2,6)*a_t - C(6,6)*a_zt;
end

% See 66_Zhao Eq 1
dBlockedTorque_dT(j) = 2*pi*trapz(r,r.^2.*stress_zt_dT(:,j).');
% See 66_Zhao Eq 12
blockedTorque(j) = ...
    trapz([T(j-1) T(j)], [dBlockedTorque_dT(j-1) dBlockedTorque_dT(j)])...
    + blockedTorque(j-1);
end

% stress_zt just for debugging and validation
% this line can be removed to improve efficiency
stress_zt = trapz(T,stress_zt_dT,2);

print = false;

end