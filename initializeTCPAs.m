function [] = initializeTCPAs(x,n)
% Initializes global tcpa nx1 structure of n_TCPAs in parallel
% This function is used to assign design variable x(i) to the appropriate
% parameter of the appropriate individual TCPA
% x  vector of all design variables for all TCPAs
% n  number of TCPAs in parallel

persistent print;
if isempty(print)
    print = true;
end

if print
    disp('> ');
    disp('> initializeTCPA.m');
    disp('> Assumed change in fiber diameter is negligible');
    print = false;
end

% Parse x to sort design variables into respective catagories
len = length(x);
if len == n % then x = [d(1) d(2) ... d(n)]
    d = x;
    C = 6*ones(n,1);
    L0 = 0.049*ones(n,1);
    L = L0 + 0.015;
    N = 46*ones(n,1);
elseif len == 2*n % then x = [d(1) ... d(n) C(1) ... C(n)];
    d = x(1:n);
    C = x(n+1:2*n);
    L0 = 0.05*ones(n,1);
    %L = L0 + 0.015;
    L = 1.2*L0;
    %N = 46*ones(n,1);
    N = L0./d;
elseif len == 4*n
    % then x = [d(1) ... d(n) C(1) ... C(n) L0(1) ... L0(n) r(1) ... r(n)]
    % where r = L/L0
    d = x(1:n);
    C = x(n+1:2*n);
    L0 = x(2*n+1:3*n);
    L = x(3*n+1:4*n).*L0;
    N = L0./d;
else
    error('x contains invalid number of design variables for n TCPAs')
end

global tcpa
tcpa = [];

T_amb = getT_amb(); % [ºC] ambient room temperature

% TCPA coil fully defined by following trigonometric relationships:
%   sin(beta_c) = L/l
%   cos(beta_c) = pi*N*D/l
% Therefore initial TCPA triangle is fully defined by L0, N, and D0; and
% stretched TCPA triangle is fully defined by L and l 

for i = 1:n
    % These values are the same for all tcpas, to easily retrieve inputs to
    % function after tcpas have been initialized
    tcpa(i).n = n;
    tcpa(i).x = x;
    
    % Independent TCPA initial parameters
    tcpa(i).d = d(i); %0.0008; % [m] fiber diam
    tcpa(i).springIndex = C(i); % defines initial coil diameter D0
    tcpa(i).N = N(i); % number of coils
    tcpa(i).L0 = L0(i); % [m] fully contracted coil length
    % beta_f_d sign defintion determines (+) homo vs (-) heterochiral
    tcpa(i).beta_f_d = 35*pi/180; % [rad] fiber bias angle from vertical @ r=d/2, CCW is +

    % Independent TCPA heating element parameters
    tcpa(i).d_wire = 0.00025; % [m] measured diam of resistance wire in lab
    tcpa(i).ironWire = false;
    % l_wire calcualted to have a resitance of 3.6 ohm @ 23ºC
    tcpa(i).l_wire = 1.8; % [m] length of resistance wire
    
    % Independent stretched TCPA parameters
    tcpa(i).L = L(i); % [m] stretched length of TCPA

    % Dependent TCPA parameters
    tcpa(i).l0 = ...
        sqrt(tcpa(i).L0^2 + pi^2*tcpa(i).N^2*tcpa(i).springIndex^2*tcpa(i).d^2);
    
    % Return error if stretched TCPA results in a coil index of D/d < 1
    % Approzimate l as l0
    D = sqrt(tcpa(i).l0^2 - tcpa(i).L^2)/(pi*tcpa(i).N); % [m] coil diameter
    if D/tcpa(i).d < 1
        error('initializeTCPAs.m: D/d < 1. TCPA is too stretched.')
    end
    
    % Calculate blocked torque, total force, and thermal force as a
    % function of temp
    [blockedTorque_vect, T_vect,~,~,delta_z] = ...
        getBlockedTorque(tcpa(i).beta_f_d,tcpa(i).d,T_amb,180);
    tcpa(i).temp = T_vect;
    tcpa(i).blockedTorque = blockedTorque_vect;
    tcpa(i).l = tcpa(i).l0*(1 + delta_z); % [m] temp-dependent fiber length
    
    for j = 1:length(T_vect)
        [F, F_th, dL_th] = getForce(tcpa(i).temp(j),tcpa(i).d,tcpa(i).l0,...
            tcpa(i).l(j),tcpa(i).L0,tcpa(i).L,tcpa(i).N,tcpa(i).blockedTorque(j));
        tcpa(i).totalForce(j) = F;
        tcpa(i).thermalForce(j) = F_th;
        tcpa(i).thermalContraction(j) = dL_th;
        
    % Calculate gain and time constant of the first-order electrothermal
    % system as a function of temperature
    	[tau, gain, R] = getElectrothermalParam(tcpa(i).temp(j),...
            tcpa(i).d,tcpa(i).l0,tcpa(i).l(j),tcpa(i).L0,tcpa(i).L,tcpa(i).N,...
            tcpa(i).d_wire,tcpa(i).l_wire,tcpa(i).ironWire);
        tcpa(i).tau(j) = tau;
        tcpa(i).gain(j) = gain;
        tcpa(i).R(j) = R;
    end
    
end

end