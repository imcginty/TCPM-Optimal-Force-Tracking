function [tau, gain, R] = getElectrothermalParam(T,d,l0,l,L0,L,N,d_wire,l_wire,ironWire)
% Returns time constant and steady state gain of first order eletrothermal
% dynamic system
% Returns NaN if stretched length stretched coil to violate D/d > 1
% T [ºC] temperature used to define air and resistive electrical properties
% d [m] fiber diameter
% l0 [m] initial length of twisted fiber
% L0 [m] initial length of TCPA 
% L [m] stretched length of TCPA 
% N number of TCPA coils
% d_wire [m] measured diam of resistance wire in lab
% l_wire [m] length of resistance wire

persistent print;
if isempty(print)
    print = true;
end
persistent printK;
if isempty(printK)
    printK = true;
end

if(print)
disp('> ')
disp('> getElectrothermalParam.m')
disp(strcat('> All properties evaluated at ',num2str(T),'ºC'))
if ironWire
    disp('> Iron wire used for Joule heating')
else
    disp('> Constantan wire used for Joule heating')
end
%disp('> Tau varies ~50% from 20-200ºC, gain varies significantly')
disp('> Electrothermal model assumes l = l0') 
end

T_amb = getT_amb(); % [ºC] ambient room temperature

if T < T_amb
error('getTemp.m: T > T_amb')
end


if T == T_amb && printK
    disp('> ')
    disp('> getElectrothermalParam.m')
	disp('> WARNING: T = T_amb. Setting T = T_amb + 1.')
    T = T_amb + 1;
    printK = false;
end

% TCPA parameters
D = sqrt(l^2 - L^2)/(pi*N); % [m] coil diameter
% D varies with stretch as L changes, and with temp as l changes
% Return NaN if coil is so stretched D/d < 1
if D/d < 1
    tau = NaN;
    gain = NaN;
    disp('> ')
    disp('> getTemp.m')
    disp(strcat('> WARNING: D/d > 1 constrain violated for d=',...
        num2str(d),'and D=',num2str(d)))
else
    
r = L/L0;

%[k Pr nu beta] = getAirProperties(T);
% Thermophysical fluid properties from "Berechnung der oberflachenbelastung
% von widerstandsheizelementen bei freier konvektion"
% Properties evaluated at mean temp Tm
Tm = 0.5*(T + T_amb); 
lambda = 2.353E-2 + 8.219E-5*Tm - 3.666E-8*Tm.^2 + 1.049E-11*Tm.^3;
nu = 1.162E-5 + 9.657E-8*Tm + 6.503E-11*Tm.^2;
Pr = 7.147E-1 - 2.452E-4*Tm + 5.819E-7*Tm.^2 - 4.702E-10*Tm.^3 + 1.382E-13*Tm.^4;
beta = 2./(T+273+getT_amb+273);

g = 9.81; % [m/s^2]

% Properties of resistance wire - from engineeringtoolbox.com
if ironWire
    resistivity = 9.71E-8; % [ohm m] resisitivity of iron at 20ºC
    resisTempCoeff = 6.41E-3; % [1/ºC] temp coeff of iron
else
    resistivity = 49E-8; % [ohm m] resisitivity of constantan at 20ºC
    resisTempCoeff = 3E-5; % [1/ºC] temp coeff of constantan
end    
T_R0 = 20; % [ºC] temperature at which resitivity is defined

% Properties of nylon - rho and c from matweb.com
rho = 1.3E6; % [g/m^3] density
c = 1.67; % [J/g ºC] specific heat capacity

% Properties of nylon - rho and c from engineeringtoolbox.com
rho = 1.14E3; % [kg/m^3] density
c = 1310; % [J/kg ºC] specific heat capacity


% Electrothermal model: C_th*dT/dt + h*A*(T-T0) = V^2/R
%   (T-T0)/V^2 = (1/(R*h*A)) / (C_th/(h*A) s + 1) 

C_th = pi/4*d^2*l0*c*rho; % [J/ºC] thermal mass of nylon, EXCLUDING IRON
% C_th = 0.357; % [Ws/ºC] for HV_C, including mass of iron

% Free convection for a coil - Hauser
% From "Berechnung der oberflachenbelastung von widerstandsheizelementen
% bei freier konvektion"
D1 = d*(1 + (pi*D/d - 1)/r);
Gr_D1 = g*beta*(T - T_amb)*D1^3/nu^2;
K = log10(Gr_D1*Pr);
if K < -6
    %disp(strcat('> WARNING: Setting K = -6. Coil free convection K < -6 for d=',...
    %    num2str(d),'; D=',num2str(D),'; T=',num2str(T)))
    K = -6;
end
if K > 4
    disp(strcat('> WARNING: Coil free convection K =',num2str(K),'> 4 for d=',...
        num2str(d),'; D=',num2str(D),'; T=',num2str(T)))
        K = 4;
end
C_D = (1.32478 + 0.03428*K + 0.00288*K^2 - 0.00003*K^4)^K;
h = C_D*lambda*r/(pi*D); % [W/m^2 K]
% h = 10.4; % [W/m^2 K] for HV_C forced convection over vertical coil
A = pi*d*l;

R = resistivity*(1 + resisTempCoeff*(T-T_R0))*l_wire/(pi/4*d_wire^2);
%R = 3.6; % [ohm] for HV_C at 23ºC

tau = C_th/(h*A); % [s] thermal time contant
gain = 1/(R*h*A); % [ºC/V^2]

end

print = false;

end
